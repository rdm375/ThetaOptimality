#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot C++ goal-specific qualification outputs; run no inference or bounds."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results" / "goal_dependence"


def read_csv(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


decisions = read_csv("decisions.csv")
envelopes = read_csv("envelopes.csv")
states = read_csv("state_reconstruction.csv")
cases = ["different_goal_times", "displacement_blocker"]
case_titles = {
    "different_goal_times": r"Different qualification times: $\theta=0.5$, $N=1024$",
    "displacement_blocker": r"One blocking goal: $\theta=0.65$, $N=4096$",
}
goals = ["displacement", "velocity", "support_force"]
goal_labels = {
    "displacement": "Displacement",
    "velocity": "Velocity",
    "support_force": "Support force",
}
colors = {"displacement": "#2876a5", "velocity": "#b45b27", "support_force": "#39824b"}

# Show the physical truth, the continuous least-squares fit to the same noisy
# samples, and the theta-discrete optimal reconstruction. The bands are the
# C++-computed nodal spreads of the shared observation-consistent set, including
# the numerical-model allowance. Each panel marks that goal's persistent time.
state_specs = [
    ("displacement", "true_displacement_m", "analytical_fit_displacement_m",
     "discrete_fit_displacement_m", "displacement_bound_m", "Displacement (m)"),
    ("velocity", "true_velocity_m_s", "analytical_fit_velocity_m_s",
     "discrete_fit_velocity_m_s", "velocity_bound_m_s", "Velocity (m/s)"),
    ("support_force", "true_support_force_N", "analytical_fit_support_force_N",
     "discrete_fit_support_force_N", "support_force_bound_N", "Support force (N)"),
]
fig, axes = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
for col, case in enumerate(cases):
    case_rows = [r for r in states if r["case"] == case]
    time = np.array([float(r["elapsed_s"]) for r in case_rows])
    case_decisions = {r["goal"]: r for r in decisions if r["case"] == case}
    for row, (goal, truth_key, analytic_key, discrete_key, bound_key, ylabel) in enumerate(state_specs):
        ax = axes[row, col]
        truth = np.array([float(r[truth_key]) for r in case_rows])
        analytic = np.array([float(r[analytic_key]) for r in case_rows])
        discrete = np.array([float(r[discrete_key]) for r in case_rows])
        bound = np.array([float(r[bound_key]) for r in case_rows])
        ax.plot(time, truth, color="black", lw=1.45, label="Analytical truth")
        ax.plot(time, analytic, color="#777777", ls="--", lw=1.1,
                label="Analytical fit to same data")
        ax.plot(time, discrete, color=colors[goal], lw=1.1,
                label="Theta-discrete fit")
        ax.fill_between(time, discrete - bound, discrete + bound,
                        color=colors[goal], alpha=.20,
                        label="Observation set + numerical allowance")
        decision = case_decisions[goal]
        if decision["qualified"] == "1":
            qtime = float(decision["persistent_time_s"])
            ax.axvline(qtime, color=colors[goal], ls=":", lw=1.35,
                       label=f"Persistent qualification: {qtime:.3f} s")
        else:
            ax.text(.98, .05, "NOT QUALIFIED", transform=ax.transAxes,
                    ha="right", va="bottom", color="#9e3226",
                    fontsize=8, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{case_titles[case]}\n{goal_labels[goal]}")
        ax.grid(alpha=.2)
        if row == 0:
            ax.legend(fontsize=7, loc="best")
        if row == len(state_specs) - 1:
            ax.set_xlabel("Elapsed time from $t_0$ (s)")
fig.suptitle("Analytical and theta-discrete reconstructions with goal-specific uncertainty and qualification")
fig.tight_layout(rect=(0, 0, 1, .95))
fig.savefig(ROOT / "goal_dependence_reconstruction.png", dpi=180)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(14, 5.3), sharey=True)
for ax, case in zip(axes, cases):
    case_decisions = {r["goal"]: r for r in decisions if r["case"] == case}
    for goal in goals:
        data = [r for r in envelopes if r["case"] == case and r["goal"] == goal
                and r["normalized_whole_interval_bound"]]
        time = np.array([float(r["elapsed_s"]) for r in data])
        ratio = np.array([float(r["normalized_whole_interval_bound"]) for r in data])
        decision = case_decisions[goal]
        status = (f"qualified from {float(decision['persistent_time_s']):.3f} s"
                  if decision["qualified"] == "1" else "not qualified")
        ax.plot(time, ratio, color=colors[goal], lw=1.2,
                label=f"{goal_labels[goal]}: {status}")
        if decision["qualified"] == "1":
            qtime = float(decision["persistent_time_s"])
            ax.axvline(qtime, color=colors[goal], ls="--", lw=.9, alpha=.8)
        else:
            ax.text(.98, .96 - .09 * goals.index(goal),
                    f"{goal_labels[goal]}: NOT QUALIFIED",
                    transform=ax.transAxes, ha="right", va="top",
                    color=colors[goal], fontsize=8,
                    bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none"})
    ax.axhline(1., color="black", ls=":", lw=1.1, label="Goal tolerance")
    joint = next(r for r in decisions if r["case"] == case)
    if joint["joint_qualified"] == "1":
        ax.axvline(float(joint["joint_time_s"]), color="black", ls="-.", lw=1.1,
                   label=f"Joint qualification: {float(joint['joint_time_s']):.3f} s")
    else:
        blocker = case_decisions["displacement"]["blocking_goals"]
        ax.text(.98, .04, f"JOINT NOT QUALIFIED — blocked by {blocker.replace('_', ' ')}",
                transform=ax.transAxes, ha="right", va="bottom", color="#9e3226",
                fontsize=8, fontweight="bold",
                bbox={"facecolor": "white", "alpha": .85, "edgecolor": "none"})
    ax.set_yscale("log")
    ax.set_xlim(0, 8)
    ax.set_xlabel("Elapsed time from $t_0$ (s)")
    ax.set_title(case_titles[case])
    ax.grid(which="both", alpha=.2)
    ax.legend(fontsize=7, loc="upper left")
axes[0].set_ylabel("Whole-interval uncertainty bound / goal tolerance")
fig.suptitle("One reconstructed state set, three engineering goals, goal-specific decisions")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "goal_dependence_bounds.png", dpi=180)
plt.close(fig)

# A compact timeline makes the conjunction rule explicit: joint qualification
# starts at the latest individual persistent time and is absent if any goal is
# unqualified.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
for ax, case in zip(axes, cases):
    case_decisions = {r["goal"]: r for r in decisions if r["case"] == case}
    for y, goal in enumerate(goals):
        decision = case_decisions[goal]
        color = colors[goal]
        ax.hlines(y, 0, 8, color="0.82", lw=4, zorder=1)
        if decision["qualified"] == "1":
            qtime = float(decision["persistent_time_s"])
            ax.hlines(y, qtime, 8, color=color, lw=5, zorder=2)
            ax.plot(qtime, y, "o", color=color, ms=6, zorder=3)
            ax.text(qtime, y + .14, f"{qtime:.3f} s", ha="center", va="bottom",
                    color=color, fontsize=8)
        else:
            ax.text(7.9, y, "NOT QUALIFIED", ha="right", va="center",
                    color="#9e3226", fontsize=8, fontweight="bold")
    joint = next(r for r in decisions if r["case"] == case)
    y = len(goals)
    ax.hlines(y, 0, 8, color="0.82", lw=4, zorder=1)
    if joint["joint_qualified"] == "1":
        qtime = float(joint["joint_time_s"])
        ax.hlines(y, qtime, 8, color="black", lw=5, zorder=2)
        ax.plot(qtime, y, "ko", ms=6, zorder=3)
        ax.text(qtime, y + .14, f"{qtime:.3f} s", ha="center", va="bottom", fontsize=8)
    else:
        blocker = joint["blocking_goals"].replace("_", " ")
        ax.text(7.9, y, f"NOT QUALIFIED ({blocker} blocks)", ha="right", va="center",
                color="#9e3226", fontsize=8, fontweight="bold")
    ax.set_yticks(range(4), [goal_labels[g] for g in goals] + ["All goals (joint)"])
    ax.invert_yaxis()
    ax.set_xlim(0, 8)
    ax.set_xlabel("Elapsed time from $t_0$ (s)")
    ax.set_title(case_titles[case])
    ax.grid(axis="x", alpha=.2)
fig.suptitle("Qualification is conjunctive: every goal must pass for the remaining horizon")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "goal_dependence_times.png", dpi=180)
plt.close(fig)
print("Plotted retained C++ goal-dependence states, uncertainty, and decisions; Python performs no inference or qualification.")
