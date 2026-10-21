#!/usr/bin/env python3
"""Plot C++ correlated-washout results; do no fitting or uncertainty analysis."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results" / "drug_washout_correlated"


def read(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


paths = read("reconstruction.csv")
budgets = read("error_budgets.csv")
decisions = read("decisions.csv")
convergence = read("convergence.csv")
comparisons = read("covariance_comparison.csv")
model_labels = {"correlated": "Correct AR(1) covariance",
                "diagonal_assumption": "Diagonal assumption (misspecified)"}
colors = {"correlated": "#2876a5", "diagonal_assumption": "#b45b27"}

# Analytical and discrete fits, with the covariance-derived confidence set and
# numerical allowance. The exact state is included to show the fixed data
# realization's truth against both inferences.
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
for model in model_labels:
    rows = [r for r in paths if r["model"] == model and int(r["steps"]) == 512]
    time = np.array([float(r["elapsed_h"]) for r in rows])
    for ax, (truth_key, analytic_key, discrete_key, bound_key, goal, ylabel) in zip(
        axes,
        (("true_concentration_mg_L", "analytical_fit_mg_L", "discrete_fit_mg_L",
          "total_concentration_bound", "concentration", "Concentration (mg/L)"),
         ("true_remaining_auc_mg_h_L", "analytical_fit_remaining_auc_mg_h_L",
          "discrete_fit_remaining_auc_mg_h_L", "total_auc_bound", "remaining_auc",
          "Remaining AUC (mg h/L)"))):
        truth = np.array([float(r[truth_key]) for r in rows])
        analytic = np.array([float(r[analytic_key]) for r in rows])
        discrete = np.array([float(r[discrete_key]) for r in rows])
        total = np.array([float(r[bound_key]) for r in rows])
        ax.plot(time, discrete, color=colors[model], lw=1.1,
                label=f"{model_labels[model]}: theta-discrete fit")
        ax.fill_between(time, discrete-total, discrete+total,
                        color=colors[model], alpha=.14,
                        label=f"{model_labels[model]}: confidence + numerical allowance")
        ax.plot(time, analytic, color=colors[model], ls="--", lw=1.0,
                label=f"{model_labels[model]}: analytical fit")
        if model == "correlated":
            ax.plot(time, truth, color="black", lw=1.5, label="Analytical truth")
        decision = next(r for r in decisions if r["model"] == model and int(r["steps"]) == 512
                        and r["goal"] == goal)
        qtime = float(decision["combined_time_h"])
        ax.axvline(qtime, color=colors[model], ls=":", lw=1.2,
                   label=f"{model_labels[model]}: qualifies {qtime:.3f} h")
        ax.set_xlabel("Elapsed time after IV bolus (h)")
        ax.set_ylabel(ylabel)
        ax.set_title("Plasma concentration" if goal == "concentration" else "Remaining exposure")
        ax.grid(alpha=.2)
axes[0].legend(fontsize=7, loc="best")
axes[1].legend(fontsize=7, loc="best")
fig.suptitle(r"Correct covariance versus an independence assumption ($N=512$, $\theta=0.5$)")
fig.tight_layout(rect=(0, 0, 1, .92))
fig.savefig(ROOT / "correlated_reconstruction.png", dpi=180)
plt.close(fig)

# Directly show how the covariance assumption changes the inferred initial
# concentration range; this avoids hiding a small parameter-space difference
# on the much larger concentration-trajectory axes.
fig, ax = plt.subplots(figsize=(8.5, 3.6))
for y, model in enumerate(("correlated", "diagonal_assumption")):
    row = next(r for r in comparisons if r["model"] == model and int(r["steps"]) == 512)
    estimate = float(row["analytical_fit_mg_L"])
    lower = float(row["confidence_lower_mg_L"])
    upper = float(row["confidence_upper_mg_L"])
    ax.errorbar(estimate, y, xerr=[[estimate-lower], [upper-estimate]], fmt="o",
                color=colors[model], capsize=5, lw=2, label=model_labels[model])
ax.axvline(10.0, color="black", ls="--", lw=1.2, label="True initial concentration")
ax.set_yticks([0, 1], [model_labels["correlated"], model_labels["diagonal_assumption"]])
ax.invert_yaxis()
ax.set_xlabel("Initial concentration $c_0$ (mg/L)")
ax.set_title("Confidence interval width changes much more than the fitted value")
ax.grid(axis="x", alpha=.25)
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .015),
           ncol=3, fontsize=8)
fig.tight_layout(rect=(0, .18, 1, 1))
fig.savefig(ROOT / "initial_concentration_intervals.png", dpi=180)
plt.close(fig)

# Compare the error components normalized by each goal tolerance. The curves
# are the already-computed C++ whole-interval bounds.
fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, sharey=True)
for row, goal in enumerate(("concentration", "remaining_auc")):
    for col, steps in enumerate((32, 512)):
        ax = axes[row, col]
        for model in model_labels:
            data = [r for r in budgets if r["model"] == model and int(r["steps"]) == steps
                    and r["goal"] == goal]
            time = np.array([float(r["elapsed_h"]) for r in data])
            tol = np.array([float(r["tolerance"]) for r in data])
            m = np.array([float(r["measurement_component"]) for r in data]) / tol
            n = np.array([float(r["numerical_component"]) for r in data]) / tol
            total = np.array([float(r["total_bound"]) for r in data]) / tol
            color = colors[model]
            ax.plot(time, m, color=color, ls="--", lw=1.0,
                    label=f"{model_labels[model]}: measurement")
            ax.plot(time, n, color=color, ls=":", lw=1.0,
                    label=f"{model_labels[model]}: numerical")
            ax.plot(time, total, color=color, lw=1.5,
                    label=f"{model_labels[model]}: total")
            decision = next(r for r in decisions if r["model"] == model and int(r["steps"]) == steps
                            and r["goal"] == goal)
            ax.axvline(float(decision["combined_time_h"]), color=color, ls="-.", lw=.9)
        ax.axhline(1., color="black", ls="--", lw=1., label="Goal tolerance")
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 5)
        ax.set_title(f"N={steps}: {'concentration' if goal == 'concentration' else 'remaining AUC'}")
        ax.grid(which="both", alpha=.2)
        if col == 0:
            ax.set_ylabel("Bound / goal tolerance")
        if row == 1:
            ax.set_xlabel("Elapsed time (h)")
axes[0, 0].legend(fontsize=6.5, loc="upper right", ncol=2)
fig.suptitle("Covariance changes the measurement uncertainty; refinement shrinks numerical error")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "correlated_error_budgets.png", dpi=180)
plt.close(fig)

# Show convergence of the complete discrete workflow against analytical
# references. Plotting displays the rate computed and written by C++.
quantity_specs = [("state_error", "State"), ("tangent_error", "Tangent"),
                  ("cost_error", "Cost"), ("gradient_error", "Gradient"),
                  ("hessian_error", "Hessian"), ("fit_error", "Fitted initial concentration"),
                  ("auc_error", "Theta-quadrature AUC")]
fig, axes = plt.subplots(2, 4, figsize=(15, 7.4))
for ax, (key, title) in zip(axes.flat, quantity_specs):
    for model, color in colors.items():
        for theta, style in ((0.5, "-"), (0.65, "--")):
            data = [r for r in convergence if r["model"] == model
                    and abs(float(r["theta"]) - theta) < 1e-12]
            h = np.array([float(r["h"]) for r in data])
            error = np.array([float(r[key]) for r in data])
            ax.loglog(h, error, marker="o", ls=style, color=color,
                      label=f"{model_labels[model]}, theta={theta:g}")
    ax.set_title(title)
    ax.set_xlabel("Step size h (h)")
    ax.set_ylabel("Absolute error")
    ax.grid(which="both", alpha=.2)
handles, labels = axes.flat[0].get_legend_handles_labels()
axes.flat[-1].axis("off")
fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=8,
           bbox_to_anchor=(.5, .005))
fig.suptitle("Correlated observation optimality-system convergence")
fig.tight_layout(rect=(0, .09, 1, .93))
fig.savefig(ROOT / "correlated_convergence.png", dpi=180)
plt.close(fig)

# Qualification-time comparison across grid refinements. These times are read
# from C++ decisions; no threshold crossing is recomputed here.
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True)
for ax, goal in zip(axes, ("concentration", "remaining_auc")):
    for model in model_labels:
        data = sorted([r for r in decisions if r["model"] == model and r["goal"] == goal],
                      key=lambda r: int(r["steps"]))
        n = np.array([int(r["steps"]) for r in data])
        t = np.array([float(r["combined_time_h"]) for r in data])
        ax.plot(n, t, "o-", color=colors[model], label=model_labels[model])
    ax.set_xscale("log", base=2)
    ax.set_xticks([32, 128, 512], ["32", "128", "512"])
    ax.set_xlabel("Theta steps N")
    ax.set_ylabel("Persistent qualification time (h)")
    ax.set_title("Concentration" if goal == "concentration" else "Remaining AUC")
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
fig.suptitle("Ignoring assay correlation makes qualification appear earlier")
fig.tight_layout(rect=(0, 0, 1, .92))
fig.savefig(ROOT / "correlated_qualification_times.png", dpi=180)
plt.close(fig)

print("Plotted C++ correlated-washout reconstruction, bounds, convergence, and decisions; Python performs no inference.")
