#!/usr/bin/env python3
"""Plot retained C++ nonuniform-mesh convergence data; run no simulation."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results" / "nonuniform_oscillator"


def read_csv(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


convergence = read_csv("convergence.csv")
mesh = read_csv("mesh.csv")
elapsed_observations = np.array([.5, 1., 2., 3., 4., 6.])
quantities = [
    ("state", "Forward state"),
    ("tangent", "Tangent linear state"),
    ("nodal_adjoint", "Nodal adjoint"),
    ("recovered_initial_state", "Recovered initial state"),
    ("running_costate", "Running-cost costate"),
]
colors = plt.cm.tab10(np.linspace(0, 1, len(quantities)))
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

selected = [row for row in mesh if np.isclose(float(row["theta"]), .5)
            and row["subdivisions"] == "16" and row["step_s"]]
step_mid = np.array([float(row["elapsed_s"]) + .5 * float(row["step_s"]) for row in selected])
step_size = np.array([float(row["step_s"]) for row in selected])
axes[0].semilogy(step_mid, step_size, color="#2876a5", lw=1.2)
for j, time in enumerate(elapsed_observations):
    axes[0].axvline(time, color="#b45b27", ls="--", lw=.8,
                    label="Fixed observation times" if j == 0 else None)
axes[0].set_xlabel("Elapsed time from $t_0$ (s)")
axes[0].set_ylabel("Step length $h_n$ (s)")
axes[0].set_title("Quadratically graded mesh, $m=16$")
axes[0].grid(which="both", alpha=.2)
axes[0].legend(fontsize=8)

for ax, theta, expected, title in [
    (axes[1], "0.5", 2, r"Midpoint: expected $O(h_{max}^2)$"),
    (axes[2], "0.65", 1, r"$θ=0.65$: expected $O(h_{max})$"),
]:
    rows = [row for row in convergence if np.isclose(float(row["theta"]), float(theta))]
    h = np.array([float(row["h_max_s"]) for row in rows])
    for color, (field, label) in zip(colors, quantities):
        error = np.array([float(row[field]) for row in rows])
        ax.loglog(h, error, "o-", color=color, label=label)
    reference_h = np.array([h.max(), h.min()])
    reference_error = .3 * (reference_h / reference_h[0]) ** expected
    ax.loglog(reference_h, reference_error, "k--", lw=1, label=f"slope {expected}")
    ax.invert_xaxis()
    ax.set_xlabel("Maximum step $h_{max}$ (s)")
    ax.set_ylabel("Scaled error")
    ax.set_title(title)
    ax.grid(which="both", alpha=.2)
    ax.legend(fontsize=7)

fig.suptitle("Theta-consistent optimality system on fixed-observation nonuniform meshes", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, .94))
fig.savefig(ROOT / "nonuniform_convergence.png", dpi=180)
plt.close(fig)

# Separate views make the mesh design, numerical quantities, and rate evidence
# readable without enlarging the overview figure.
fig, ax = plt.subplots(figsize=(9, 4.8))
ax.semilogy(step_mid, step_size, color="#2876a5", lw=1.2, label="Step length")
for j, time in enumerate(elapsed_observations):
    ax.axvline(time, color="#b45b27", ls="--", lw=.8,
               label="Fixed observation times" if j == 0 else None)
ax.set_xlabel("Elapsed time from $t_0$ (s)")
ax.set_ylabel("Step length $h_n$ (s)")
ax.set_title("Quadratically graded mesh, $m=16$: observation nodes stay fixed")
ax.grid(which="both", alpha=.2)
ax.legend()
fig.tight_layout()
fig.savefig(ROOT / "mesh_spacing.png", dpi=180)
plt.close(fig)

state_quantities = [
    ("state", "Forward state"),
    ("tangent", "Tangent linear state"),
    ("nodal_adjoint", "Nodal adjoint"),
    ("recovered_initial_state", "Recovered initial state"),
]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
for ax, theta, expected in zip(axes, ("0.5", "0.65"), (2, 1)):
    rows = [row for row in convergence if np.isclose(float(row["theta"]), float(theta))]
    h = np.array([float(row["h_max_s"]) for row in rows])
    for color, (field, label) in zip(plt.cm.tab10(np.linspace(0, 1, len(state_quantities))), state_quantities):
        ax.loglog(h, [float(row[field]) for row in rows], "o-", color=color, label=label)
    ax.invert_xaxis()
    ax.set_xlabel("Maximum step $h_{max}$ (s)")
    ax.set_ylabel("Scaled error")
    ax.set_title(f"$\\theta={theta}$; expected order {expected}")
    ax.grid(which="both", alpha=.2)
    ax.legend(fontsize=8)
fig.suptitle("State, tangent, adjoint, and recovered-state convergence")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "state_tangent_adjoint_convergence.png", dpi=180)
plt.close(fig)

objective_quantities = [
    ("point_cost", "Point-data cost"),
    ("point_gradient", "Point-data gradient"),
    ("point_hessian", "Point-data Hessian"),
    ("running_cost", "Running cost"),
    ("running_gradient", "Running gradient"),
    ("running_costate", "Running costate"),
    ("running_hessian", "Running Hessian"),
]
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
for ax, theta, expected in zip(axes, ("0.5", "0.65"), (2, 1)):
    rows = [row for row in convergence if np.isclose(float(row["theta"]), float(theta))]
    h = np.array([float(row["h_max_s"]) for row in rows])
    for color, (field, label) in zip(plt.cm.tab10(np.linspace(0, 1, len(objective_quantities))), objective_quantities):
        ax.loglog(h, [float(row[field]) for row in rows], "o-", color=color, label=label)
    ax.invert_xaxis()
    ax.set_xlabel("Maximum step $h_{max}$ (s)")
    ax.set_ylabel("Relative error")
    ax.set_title(f"$\\theta={theta}$; expected order {expected}")
    ax.grid(which="both", alpha=.2)
    ax.legend(fontsize=7)
fig.suptitle("Point-observation and running-cost optimality quantities")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "objective_derivative_convergence.png", dpi=180)
plt.close(fig)

order_rows = read_csv("orders.csv")
quantity_order = [row[0] for row in quantities] + [
    "point_cost", "point_gradient", "point_hessian", "running_cost",
    "running_gradient", "running_hessian",
]
labels = [
    "State", "Tangent", "Nodal adjoint", "Recovered state", "Running costate",
    "Point cost", "Point gradient", "Point Hessian", "Running cost",
    "Running gradient", "Running Hessian",
]
fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
for ax, theta, expected in zip(axes, ("0.5", "0.65"), (2, 1)):
    rows = {row["quantity"]: row for row in order_rows
            if np.isclose(float(row["theta"]), float(theta))}
    values = [float(rows[key]["observed_order"]) for key in quantity_order]
    passed = [rows[key]["passed"] == "1" for key in quantity_order]
    colors = ["#2876a5" if result else "#b45b27" for result in passed]
    ax.scatter(np.arange(len(labels)), values, color=colors, zorder=3)
    ax.axhline(expected, color="black", ls="--", lw=1, label=f"Expected order {expected}")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=55, ha="right")
    ax.set_ylim(.75, 2.2)
    ax.set_ylabel("Observed order from final mesh pair")
    ax.set_title(f"$\\theta={theta}$")
    ax.grid(axis="y", alpha=.2)
    ax.legend()
fig.suptitle("Measured rates across all eleven quantities (blue = within test tolerance)")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "observed_orders.png", dpi=180)
plt.close(fig)

# State reconstruction with the deterministic physical-set uncertainty envelope.
# The bands include both the observation-consistent initial-state set and the
# propagated numerical-model allowance calculated by the C++ experiment.
state_rows = read_csv("state_uncertainty.csv")
qualification_rows = read_csv("qualification.csv")
state_specs = [
    ("true_displacement_m", "center_displacement_m", "displacement_bound_m",
     "Displacement (m)", "displacement"),
    ("true_velocity_m_s", "center_velocity_m_s", "velocity_bound_m_s",
     "Velocity (m/s)", "velocity"),
    ("true_support_force_N", "center_support_force_N", "support_force_bound_N",
     "Support force (N)", "support_force"),
]
fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for row_index, theta in enumerate(("0.5", "0.65")):
    data = [r for r in state_rows if np.isclose(float(r["theta"]), float(theta))]
    time = np.array([float(r["elapsed_s"]) for r in data])
    for col_index, (truth_key, center_key, bound_key, ylabel, goal) in enumerate(state_specs):
        ax = axes[row_index, col_index]
        truth = np.array([float(r[truth_key]) for r in data])
        center = np.array([float(r[center_key]) for r in data])
        bound = np.array([float(r[bound_key]) for r in data])
        ax.plot(time, truth, color="black", lw=1.5, label="Analytical truth")
        ax.plot(time, center, color="#2876a5", lw=1.1, label="Discrete reconstruction")
        ax.fill_between(time, center - bound, center + bound, color="#69a9cf",
                        alpha=.28, label="Qualified-set + numerical uncertainty")
        q = next(r for r in qualification_rows
                 if np.isclose(float(r["theta"]), float(theta)) and r["goal"] == goal)
        if q["qualified"] == "1":
            qt = float(q["persistent_time_s"])
            ax.axvline(qt, color="#b45b27", ls="--", lw=1.2,
                       label=f"Persistent qualification: {qt:g} s")
        else:
            ax.text(.98, .04, "NOT QUALIFIED", transform=ax.transAxes, ha="right",
                    va="bottom", color="#a23a25", fontsize=9, fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.set_title(f"$\\theta={theta}$ — {goal.replace('_', ' ')}")
        ax.grid(alpha=.22)
        if row_index == 0 and col_index == 0:
            ax.legend(fontsize=7, loc="best")
        if row_index == 1:
            ax.set_xlabel("Elapsed time from $t_0$ (s)")
fig.suptitle("Nonuniform-mesh state reconstruction, uncertainty, and persistent qualification")
fig.tight_layout(rect=(0, 0, 1, .95))
fig.savefig(ROOT / "state_uncertainty_qualification.png", dpi=180)
plt.close(fig)

# Whole-interval bounds are the actual qualification criterion; node bounds are
# shown for context. A ratio below one means that the corresponding tolerance is
# established. The horizontal qualification line marks the start of a persistent
# suffix, not merely the first isolated time at which a bound dips below one.
envelope_rows = read_csv("qualification_envelopes.csv")
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharex=True)
for ax, (_, _, _, ylabel, goal) in zip(axes, state_specs):
    for theta, color, style in (("0.5", "#2876a5", "-"), ("0.65", "#b45b27", "-")):
        data = [r for r in envelope_rows if np.isclose(float(r["theta"]), float(theta))
                and r["goal"] == goal and r["whole_interval_bound"]]
        time = np.array([float(r["elapsed_s"]) for r in data])
        tolerance = float(data[0]["tolerance"])
        node_ratio = np.array([float(r["node_bound"]) / tolerance for r in data])
        interval_ratio = np.array([float(r["whole_interval_bound"]) / tolerance for r in data])
        ax.plot(time, interval_ratio, color=color, ls=style, lw=1.25,
                label=fr"$\theta={theta}$ whole interval")
        ax.plot(time, node_ratio, color=color, ls=":", lw=1,
                label=fr"$\theta={theta}$ nodes")
        q = next(r for r in qualification_rows
                 if np.isclose(float(r["theta"]), float(theta)) and r["goal"] == goal)
        if q["qualified"] == "1":
            ax.axvline(float(q["persistent_time_s"]), color=color, ls="--", lw=.9)
    ax.axhline(1., color="black", ls="--", lw=1, label="Tolerance" if goal == "displacement" else None)
    ax.set_title(goal.replace("_", " "))
    ax.set_xlabel("Elapsed time from $t_0$ (s)")
    ax.set_ylabel("Uncertainty bound / goal tolerance")
    ax.set_ylim(bottom=0)
    ax.grid(alpha=.22)
axes[0].legend(fontsize=7, loc="best")
fig.suptitle("Persistent qualification uses the whole-interval uncertainty bound")
fig.tight_layout(rect=(0, 0, 1, .93))
fig.savefig(ROOT / "qualification_bounds.png", dpi=180)
plt.close(fig)
print("Plotted retained C++ nonuniform-mesh outputs; Python performs no simulation or order estimation.")
