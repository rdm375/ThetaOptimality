#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot C++ bounded-washout results; this script performs no inference."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "bounded_washout"

def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

summary = rows(OUT / "summary.csv")
curves = rows(OUT / "trajectory_qualification.csv")
by_case = {}
for row in curves:
    by_case.setdefault(row["case"], []).append(row)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
for ax, info in zip(axes, summary):
    name = info["case"]
    data = by_case[name]
    t = [float(r["elapsed_h"]) for r in data]
    truth = [float(r["truth_mg_L"]) for r in data]
    discrete = [float(r["discrete_fit_mg_L"]) for r in data]
    center = [float(r["feasible_center_mg_L"]) for r in data]
    width = [float(r["uncertainty_bound_mg_L"]) for r in data]
    ax.fill_between(t, [x-w for x, w in zip(center, width)],
                    [x+w for x, w in zip(center, width)], alpha=.22,
                    color="tab:blue", label="bounded-data + numerical allowance")
    ax.plot(t, truth, "k--", label="synthetic truth")
    ax.plot(t, discrete, color="tab:orange", label="theta-discrete constrained fit")
    ax.axhline(0., color="0.5", linewidth=.8)
    qtime = float(info["qualification_time_h"])
    if qtime >= 0:
        ax.axvline(qtime, color="tab:green", linestyle=":", label=f"qualified from {qtime:.2f} h")
    active = info["active_bound"]
    ax.set_title(f"{name.replace('_', ' ')}\nactive bound: {active}")
    ax.set_xlabel("elapsed time (h)")
    ax.grid(True, alpha=.25)
axes[0].set_ylabel("concentration (mg/L)")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, -.08), ncol=4)
fig.suptitle("Bounded initial concentration and whole-interval qualification")
fig.savefig(OUT / "state_qualification.png", dpi=180, bbox_inches="tight")
plt.close(fig)

refinement = rows(OUT / "refinement.csv")
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
for name in by_case:
    data = [r for r in refinement if r["case"] == name]
    h = [float(r["dt_h"]) for r in data]
    state = [max(float(r["maximum_state_error"]), 1e-18) for r in data]
    if any(float(r["maximum_state_error"]) > 1e-16 for r in data):
        axes[0].loglog(h, state, "o-", label=name.replace("_", " "))
interior = [r for r in refinement if r["case"] == "interior"]
h = [float(r["dt_h"]) for r in interior]
state_ref = [max(float(r["maximum_state_error"]), 1e-18) for r in interior]
control_ref = [max(float(r["control_error"]), 1e-18) for r in interior]
tangent_ref = [float(r["maximum_tangent_error"]) for r in interior]
axes[0].loglog(h, [state_ref[0]*(x/h[0])**2 for x in h], "k--", label="$O(h^2)$ reference")
axes[1].loglog(h, tangent_ref, "o-", label="discrete tangent")
axes[2].loglog(h, control_ref, "o-", label="interior fitted initial concentration")
axes[2].loglog(h, [control_ref[0]*(x/h[0])**2 for x in h], "k--", label="$O(h^2)$ reference")
axes[0].text(.04, .06, "Lower-bound state is exactly zero", transform=axes[0].transAxes)
axes[0].set(title="Forward state error", xlabel="step size h (h)", ylabel="maximum absolute error (mg/L)")
axes[1].set(title="Tangent error", xlabel="step size h (h)", ylabel="maximum absolute sensitivity error")
axes[2].set(title="Recovered interior control error", xlabel="step size h (h)", ylabel="absolute error (mg/L)")
for ax in axes:
    ax.grid(True, which="both", alpha=.25)
    ax.legend()
fig.suptitle("Midpoint convergence against exact washout references")
fig.savefig(OUT / "theta_convergence.png", dpi=180)
plt.close(fig)

names = [r["case"].replace("_", " ") for r in summary]
lower = [float(r["lower_multiplier"]) for r in summary]
upper = [float(r["upper_multiplier"]) for r in summary]
fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
index = list(range(len(names)))
ax.bar([i-.18 for i in index], lower, width=.36, label="$\\lambda_L$")
ax.bar([i+.18 for i in index], upper, width=.36, label="$\\lambda_U$")
ax.set_xticks(index, names, rotation=12)
ax.set_ylabel("KKT multiplier (objective units per normalized control)")
ax.set_title("Active-bound multipliers distinguish constrained optima")
ax.grid(True, axis="y", alpha=.25)
ax.legend()
fig.savefig(OUT / "kkt_multipliers.png", dpi=180)
plt.close(fig)
print("Plotted bounded-washout C++ outputs; Python performs no optimization or qualification.")
