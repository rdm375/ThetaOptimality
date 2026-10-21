#!/usr/bin/env python3
"""Render C++ drug-washout reconstruction, budget, and convergence outputs.

This script performs plotting only. It does not fit parameters, estimate error
bounds, calculate convergence rates, or make qualification decisions.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results" / "drug_washout"


def read_csv(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


reconstruction = read_csv("reconstruction.csv")
budgets = read_csv("error_budgets.csv")
decisions = read_csv("decisions.csv")
convergence = read_csv("convergence.csv")
orders = read_csv("orders.csv")

case_specs = [
    ("N32_precise", r"$N=32$, sensor bound $\pm0.005$ mg/L"),
    ("N512_precise", r"$N=512$, sensor bound $\pm0.005$ mg/L"),
    ("N512_coarse_sensors", r"$N=512$, sensor bound $\pm0.050$ mg/L"),
]
goal_specs = {
    "concentration": (
        "true_concentration_mg_L", "analytical_fit_mg_L", "discrete_fit_mg_L",
        "measurement_bound_mg_L", "numerical_bound_mg_L", "total_bound_mg_L",
        "true_concentration_mg_L", "mg/L"),
    "remaining_auc": (
        "true_remaining_auc_mg_h_L", "analytical_fit_remaining_auc_mg_h_L",
        "discrete_fit_remaining_auc_mg_h_L", "measurement_auc_bound_mg_h_L",
        "numerical_auc_bound_mg_h_L", "total_auc_bound_mg_h_L",
        "true_remaining_auc_mg_h_L", "mg h/L"),
}
colors = {"concentration": "#2876a5", "remaining_auc": "#b45b27"}
decision_by_case = {r["case"]: r for r in decisions}

# State and remaining-exposure reconstructions: physical truth and both fitted
# paths are shown with the additive sensor/model-error allowance. The lines
# and bands at nodes are connected linearly for display and whole-interval
# qualification, as documented in the experiment.
fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for col, (case, title) in enumerate(case_specs):
    rows = [r for r in reconstruction if r["case"] == case]
    time = np.array([float(r["elapsed_h"]) for r in rows])
    decision = decision_by_case[case]
    for row, goal in enumerate(("concentration", "remaining_auc")):
        ax = axes[row, col]
        (truth_key, analytic_key, discrete_key, measurement_key, numerical_key,
         total_key, _, unit) = goal_specs[goal]
        truth = np.array([float(r[truth_key]) for r in rows])
        analytic = np.array([float(r[analytic_key]) for r in rows])
        discrete = np.array([float(r[discrete_key]) for r in rows])
        bound = np.array([float(r[total_key]) for r in rows])
        color = colors[goal]
        ax.plot(time, truth, color="black", lw=1.5, label="Analytical truth")
        ax.plot(time, analytic, color="#777777", ls="--", lw=1.1,
                label="Analytical fit to data")
        ax.plot(time, discrete, color=color, lw=1.15, label="Theta-discrete fit")
        ax.fill_between(time, discrete - bound, discrete + bound, color=color,
                        alpha=.20, label="Sensor + numerical allowance")
        if goal == "concentration":
            qtime = float(decision["combined_time_h"])
            tolerance = float(decision["concentration_tolerance_mg_L"])
        else:
            qtime = float(decision["combined_auc_time_h"])
            tolerance = float(decision["remaining_auc_tolerance_mg_h_L"])
        ax.axvline(qtime, color=color, ls=":", lw=1.25,
                   label=f"Persistent qualification: {qtime:.3g} h")
        ax.set_title(f"{title}\n{'Concentration' if goal == 'concentration' else 'Remaining AUC'}")
        ax.set_ylabel(f"{'Concentration' if goal == 'concentration' else 'Remaining AUC'} ({unit})")
        ax.grid(alpha=.2)
        if col == 0 and row == 0:
            ax.legend(fontsize=7, loc="best")
        if row == 1:
            ax.set_xlabel("Elapsed time after IV bolus (h)")
        # Show the engineering tolerance for scale without implying the
        # tolerance is an uncertainty band.
        ax.text(.98, .06, f"goal tolerance = {tolerance:g} {unit}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
fig.suptitle("Drug washout: analytical and theta-discrete reconstructions with error allowances")
fig.tight_layout(rect=(0, 0, 1, .94))
fig.savefig(ROOT / "reconstruction_and_qualification.png", dpi=180)
plt.close(fig)

# Budget decomposition in units of each goal tolerance. All bounds and
# qualification times were calculated and written by the C++ executable.
fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True, sharey=True)
for col, (case, title) in enumerate(case_specs):
    case_rows = [r for r in budgets if r["case"] == case]
    decision = decision_by_case[case]
    for row, goal in enumerate(("concentration", "remaining_auc")):
        ax = axes[row, col]
        data = [r for r in case_rows if r["goal"] == goal]
        time = np.array([float(r["elapsed_h"]) for r in data])
        tolerance = np.array([float(r["tolerance"]) for r in data])
        meas = np.array([float(r["measurement_component"]) for r in data]) / tolerance
        numerical = np.array([float(r["numerical_component"]) for r in data]) / tolerance
        total = np.array([float(r["total_bound"]) for r in data]) / tolerance
        ax.plot(time, meas, color="#39824b", lw=1.2, label="Measurement component")
        ax.plot(time, numerical, color="#b45b27", lw=1.2, label="Numerical component")
        ax.plot(time, total, color="#2876a5", lw=1.5, label="Sum used for decision")
        ax.axhline(1., color="black", ls=":", lw=1.1, label="Goal tolerance")
        qtime = float(decision["combined_time_h"] if goal == "concentration"
                      else decision["combined_auc_time_h"])
        ax.axvline(qtime, color="#2876a5", ls="--", lw=1.0)
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 20)
        ax.set_title(f"{title}\n{'Concentration' if goal == 'concentration' else 'Remaining AUC'}")
        ax.grid(which="both", alpha=.2)
        if col == 0:
            ax.set_ylabel("Error allowance / goal tolerance")
        if row == 1:
            ax.set_xlabel("Elapsed time (h)")
        if col == 0 and row == 0:
            ax.legend(fontsize=7, loc="upper right")
fig.suptitle("Measurement and time-discretization error contributions")
fig.tight_layout(rect=(0, 0, 1, .94))
fig.savefig(ROOT / "error_budget_decomposition.png", dpi=180)
plt.close(fig)

# Convergence of the discrete forward/optimality system against analytical
# references. The rates below come from C++ orders.csv; the renderer only
# displays them and the precomputed error data.
quantities = [
    ("state_error", "State"), ("tangent_error", "Tangent"),
    ("cost_error", "Cost"), ("gradient_error", "Gradient"),
    ("hessian_error", "Hessian"), ("recovered_initial_error", "Recovered initial concentration"),
    ("auc_error", "Theta-quadrature AUC"),
]
fig, axes = plt.subplots(2, 4, figsize=(15, 7.2))
for ax, (column, title) in zip(axes.flat, quantities):
    for theta, color, label in (("0.5", "#2876a5", r"$\theta=0.5$"),
                                ("0.65", "#b45b27", r"$\theta=0.65$")):
        rows = [r for r in convergence if abs(float(r["theta"]) - float(theta)) < 1e-12]
        h = np.array([float(r["h"]) for r in rows])
        err = np.array([float(r[column]) for r in rows])
        rate = next(float(r["observed_order"]) for r in orders
                    if abs(float(r["theta"]) - float(theta)) < 1e-12 and
                    r["quantity"] == column.removesuffix("_error"))
        ax.loglog(h, err, "o-", color=color,
                  label=f"{label}; final order {rate:.3f}")
    ax.set_title(title)
    ax.set_xlabel("Step size $h$ (h)")
    ax.set_ylabel("Absolute error")
    ax.grid(which="both", alpha=.2)
    ax.legend(fontsize=7)
axes.flat[-1].axis("off")
fig.suptitle("Theta-system convergence against exact washout references")
fig.tight_layout(rect=(0, 0, 1, .94))
fig.savefig(ROOT / "convergence.png", dpi=180)
plt.close(fig)

# Persistent qualification timeline separates the two goals and shows how
# refining the numerical model differs from tightening the measurements.
fig, axes = plt.subplots(figsize=(12, 5.7))
for y, (case, title) in enumerate((
    ("N32_precise", "N=32, precise sensors"),
    ("N128_precise", "N=128, precise sensors"),
    ("N512_precise", "N=512, precise sensors"),
    ("N32_coarse_sensors", "N=32, coarse sensors"),
    ("N128_coarse_sensors", "N=128, coarse sensors"),
    ("N512_coarse_sensors", "N=512, coarse sensors"),
)):
    decision = decision_by_case[case]
    for offset, (goal, key, color) in enumerate((
        ("Concentration", "combined_time_h", colors["concentration"]),
        ("Remaining AUC", "combined_auc_time_h", colors["remaining_auc"]),
    )):
        qtime = float(decision[key])
        yy = y + (offset - .5) * .28
        axes.hlines(yy, 0, 8, color="0.86", lw=3)
        axes.hlines(yy, qtime, 8, color=color, lw=4)
        axes.plot(qtime, yy, "o", color=color, ms=5)
        axes.text(qtime, yy + .055, f"{qtime:.3g}", ha="center", va="bottom",
                  fontsize=7, color=color)
axes.set_yticks(range(6), ["N=32, precise sensors", "N=128, precise sensors",
                           "N=512, precise sensors", "N=32, coarse sensors",
                           "N=128, coarse sensors", "N=512, coarse sensors"])
axes.invert_yaxis()
axes.set_xlim(0, 8)
axes.set_xlabel("Elapsed time after IV bolus (h)")
axes.set_title("Persistent qualification times (concentration blue; remaining AUC orange)")
axes.grid(axis="x", alpha=.2)
fig.tight_layout()
fig.savefig(ROOT / "qualification_times.png", dpi=180)
plt.close(fig)

print("Plotted C++ drug-washout reconstructions, error budgets, convergence, and qualification times; Python performs no inference or qualification.")
