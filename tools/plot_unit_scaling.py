#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot C++ unit-conversion checks; this script performs no inference."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results" / "unit_scaling"


def read(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


trajectory = read("trajectory_comparison.csv")
comparison = read("unit_comparison.csv")
orders = read("orders.csv")
qualification = read("qualification_comparison.csv")

# Converted physical trajectories lie on top of each other. The lower panel
# shows their difference, which exposes implementation and roundoff errors.
time = np.array([float(r["time_h"]) for r in trajectory])
base = np.array([float(r["mg_L"]) for r in trajectory])
converted = np.array([float(r["ug_L_converted_to_mg_L"]) for r in trajectory])
difference = np.array([float(r["absolute_difference_mg_L"]) for r in trajectory])
fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                         gridspec_kw={"height_ratios": [3, 1]})
axes[0].plot(time, base, color="#2876a5", lw=2,
             label="mg/L, hours")
axes[0].plot(time, converted, color="#b45b27", ls="--", lw=1.4,
             label="microgram/L, minutes converted back")
axes[0].set_ylabel("Concentration (mg/L)")
axes[0].set_title("The converted discrete state is the same physical trajectory")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=.25)
axes[1].semilogy(time, np.maximum(difference, 1e-18), color="#5a7d35")
axes[1].set_xlabel("Elapsed time after IV bolus (h)")
axes[1].set_ylabel("Absolute difference (mg/L)")
axes[1].grid(which="both", alpha=.25)
fig.tight_layout()
fig.savefig(ROOT / "trajectory_unit_conversion.png", dpi=180)
plt.close(fig)

# Compare the transformed quantities across all three meshes. Values are
# relative differences after converting back to the mg/L-hour representation.
names = ["fit_initial_concentration", "objective_at_diagnostic_state",
         "scaled_initial_gradient", "scaled_initial_hessian", "total_auc",
         "auc_initial_sensitivity"]
labels = ["Fitted $c_0$", "Objective", "Scaled gradient", "Scaled Hessian",
          "Total AUC", "AUC sensitivity"]
fig, ax = plt.subplots(figsize=(10, 5.5))
colors = {32: "#b45b27", 128: "#2876a5", 512: "#5a7d35"}
for offset, quantity in enumerate(names):
    rows = [r for r in comparison if r["quantity"] == quantity]
    for row in rows:
        steps = int(row["steps"])
        error = max(float(row["relative_difference"]), 1e-18)
        ax.scatter(offset + (steps - 128) * .00035, error, s=42,
                   color=colors[steps], label=f"N={steps}" if offset == 0 else None)
ax.set_yscale("log")
ax.set_xticks(range(len(names)), labels, rotation=15, ha="right")
ax.set_ylabel("Relative difference after unit conversion")
ax.set_title("Objective and derivative quantities obey their unit transformations")
ax.grid(axis="y", which="both", alpha=.25)
ax.legend(fontsize=8, ncol=3)
fig.tight_layout()
fig.savefig(ROOT / "unit_transformed_quantities.png", dpi=180)
plt.close(fig)

# Both systems show the same order after errors are expressed in the same
# physical units. The small horizontal offsets keep the markers distinguishable.
quantity_order = ["state", "tangent", "cost", "gradient", "hessian", "fit", "auc"]
display = ["State", "Tangent", "Cost", "Gradient", "Hessian", "Fit", "AUC"]
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
unit_specs = [("mg_L_hour", "mg/L and hours", -.08, "#2876a5"),
              ("ug_L_min", "microgram/L and minutes", .08, "#b45b27")]
for ax, theta in zip(axes, (0.5, 0.65)):
    for unit, label, dx, color in unit_specs:
        rows = [r for r in orders if r["units"] == unit and float(r["theta"]) == theta]
        values = [float(next(r["observed_order"] for r in rows if r["quantity"] == q))
                  for q in quantity_order]
        ax.scatter(np.arange(len(values)) + dx, values, color=color, s=36, label=label)
    expected = 2.0 if theta == .5 else 1.0
    ax.axhline(expected, color="black", ls="--", lw=1,
               label=f"Expected order {expected:g}")
    ax.set_xticks(range(len(display)), display, rotation=25, ha="right")
    ax.set_title(rf"$\theta={theta:g}$")
    ax.grid(axis="y", alpha=.25)
    ax.legend(fontsize=7)
axes[0].set_ylabel("Observed order, last mesh pair")
fig.suptitle("Unit conversion preserves temporal convergence")
fig.tight_layout()
fig.savefig(ROOT / "unit_scaling_convergence.png", dpi=180)
plt.close(fig)

# Persistent qualification times are plotted in hours after the minute-based
# results have been converted back. Overlapping markers are the desired result.
fig, ax = plt.subplots(figsize=(8.5, 4.8))
goal_style = {"concentration": ("Concentration", "o", "#2876a5"),
              "remaining_auc": ("Remaining AUC", "s", "#b45b27")}
for goal, (label, marker, color) in goal_style.items():
    rows = [r for r in qualification if r["goal"] == goal]
    mesh = np.array([int(r["steps"]) for r in rows])
    base_time = np.array([float(r["mg_L_hour_time_h"]) for r in rows])
    converted_time = np.array([float(r["converted_time_h"]) for r in rows])
    ax.plot(mesh, base_time, marker=marker, color=color, lw=1.8,
            label=f"{label}: mg/L-hours")
    ax.scatter(mesh, converted_time, marker=marker, facecolors="none",
               edgecolors=color, s=75, label=f"{label}: converted microgram/L-minutes")
ax.set_xscale("log", base=2)
ax.set_xticks([32, 128, 512], ["32", "128", "512"])
ax.set_xlabel("Theta-method steps $N$")
ax.set_ylabel("Persistent qualification time (h)")
ax.set_title("Qualification decisions are unchanged by unit choice")
ax.grid(alpha=.25)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(ROOT / "unit_scaling_qualification.png", dpi=180)
plt.close(fig)

print("Plotted C++ unit-scaling outputs; Python performs no modeling, fitting, or qualification.")
