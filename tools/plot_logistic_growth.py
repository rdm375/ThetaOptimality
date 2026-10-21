#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot retained C++ logistic-growth state-finding and optimality outputs."""
from pathlib import Path
import csv

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "logistic_growth"


def read_rows(name):
    with (OUT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


summary = read_rows("summary.csv")
curves = read_rows("qualification_curves.csv")
intervals = read_rows("qualification_intervals.csv")
observations = read_rows("observations.csv")
fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
for col, case in enumerate(summary):
    ax = axes[0, col]
    qax = axes[1, col]
    n = int(case["steps"])
    data = [row for row in curves if int(row["steps"]) == n]
    time = [float(row["elapsed_h"]) for row in data]
    truth = [float(row["truth_g_L"]) for row in data]
    discrete = [float(row["theta_fit_g_L"]) for row in data]
    analytic_fit = [float(row["analytic_fit_g_L"]) for row in data]
    lower = [float(row["feasible_lower_g_L"]) for row in data]
    upper = [float(row["feasible_upper_g_L"]) for row in data]
    ax.fill_between(time, lower, upper, color="tab:blue", alpha=.2,
                    label="bounded-data set + theta allowance")
    ax.plot(time, truth, "k--", label="synthetic truth")
    ax.plot(time, discrete, color="tab:orange", label="theta-discrete fit")
    ax.plot(time, analytic_fit, color="tab:green", linestyle=":",
            label="analytical least-squares trajectory")
    ax.errorbar([float(row["elapsed_h"]) for row in observations],
                [float(row["measured_biomass_g_L"]) for row in observations],
                yerr=[float(row["assay_bound_g_L"]) for row in observations],
                fmt="ko", markersize=3, capsize=2, label="assays with bounded error")
    qtime = (float(case["qualification_time_h"])
             if case["qualification_time_h"] != "unqualified" else None)
    if qtime is not None:
        ax.axvline(qtime, color="tab:red", linestyle="--",
                   label=f"persistent qualification: {qtime:.2f} h")
    status = f"qualified from {qtime:.2f} h" if qtime is not None else "not qualified"
    ax.set(title=f"State reconstruction, N={n} ({status})",
           ylabel="biomass concentration (g/L)")
    ax.grid(True, alpha=.25)
    interval_data = [row for row in intervals if int(row["steps"]) == n]
    qtime_series = [float(row["interval_start_h"]) for row in interval_data]
    qbound = [float(row["total_bound_g_L"]) for row in interval_data]
    qax.plot(qtime_series, qbound, color="tab:purple", label="whole-step error bound")
    qax.axhline(.02, color="tab:red", linestyle="--", label="goal tolerance")
    if qtime is not None:
        qax.axvline(qtime, color="tab:green", linestyle=":",
                    label=f"qualified from {qtime:.2f} h")
    qax.set(title="Persistent qualification bound", xlabel="elapsed time (h)",
            ylabel="error bound (g/L)", xlim=(0, 4))
    qax.grid(True, alpha=.25)
axes[0, 0].legend(loc="upper left", fontsize=8)
axes[1, 0].legend(loc="upper right", fontsize=8)
axes[1, 1].legend(loc="upper right", fontsize=8)
for ax in axes[0]:
    ax.set_xlabel("elapsed time (h)")
fig.suptitle("Logistic-growth reconstruction and persistent state qualification")
fig.savefig(OUT / "state_qualification.png", dpi=180, bbox_inches="tight")
plt.close(fig)

convergence = read_rows("convergence.csv")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
panels = [
    (axes[0], [("state", "state"), ("tangent", "tangent"), ("adjoint", "costate")],
     "State and sensitivity convergence"),
    (axes[1], [("cost", "objective"), ("gradient", "gradient"), ("hessian", "Hessian")],
     "Discrete optimality convergence"),
    (axes[2], [("recovered_initial", "recovered initial state")],
     "Inverse-solution convergence"),
]
for ax, quantities, title in panels:
    for theta in (0.5, 0.65):
        data = sorted((row for row in convergence if float(row["theta"]) == theta),
                      key=lambda row: float(row["h"]))
        h = [float(row["h"]) for row in data]
        for column, label in quantities:
            error = [max(float(row[column]), 1e-18) for row in data]
            ax.loglog(h, error, "o-", label=f"$\\theta={theta}$ {label}")
            order = 2 if theta == .5 else 1
            ax.loglog(h, [error[0]*(step/h[0])**order for step in h], "k:", alpha=.3)
    ax.set(title=title, xlabel="step size h (h)", ylabel="absolute error")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=8)
fig.suptitle("Observed convergence of the nonlinear discrete optimality system")
fig.savefig(OUT / "theta_convergence.png", dpi=180)
plt.close(fig)

optimality = read_rows("optimality_curves.csv")
u = [float(row["control_u"]) for row in optimality]
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
specs = [
    ("cost", "exact_cost", "cost_N32", "cost_N128", "objective value"),
    ("gradient", "exact_gradient", "gradient_N32", "gradient_N128", "objective gradient"),
    ("hessian", "exact_hessian", "hessian_N32", "hessian_N128", "objective Hessian"),
]
for ax, (name, exact, coarse, fine, ylabel) in zip(axes, specs):
    for key, label in ((exact, "analytical"), (coarse, "$N=32$"), (fine, "$N=128$")):
        ax.plot(u, [float(row[key]) for row in optimality], label=label)
    ax.axvline(float(summary[1]["analytical_fit_g_L"])/8., color="0.35", linestyle=":",
               label="analytical minimizer")
    ax.set(xlabel="normalized initial biomass $u_0$", ylabel=ylabel)
    ax.grid(True, alpha=.25)
    ax.legend()
fig.suptitle("Analytical and theta-discrete objective derivatives")
fig.savefig(OUT / "objective_derivatives.png", dpi=180)
plt.close(fig)
print("Plotted retained logistic-growth C++ outputs; Python performs no simulation or inference.")
