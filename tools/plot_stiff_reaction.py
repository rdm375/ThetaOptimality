#!/usr/bin/env python3
"""Plot retained C++ results for the analytically solvable stiff reaction case."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", default="results/stiff_reaction", type=Path)
    args = parser.parse_args()
    root = args.results
    diagnostic = read_csv(root / "stiffness_trajectories.csv")
    trajectories = read_csv(root / "trajectories.csv")
    observations = read_csv(root / "observations.csv")
    convergence = read_csv(root / "convergence.csv")
    orders = read_csv(root / "orders.csv")
    qualification = read_csv(root / "qualification.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))
    for theta in (0.5, 0.65, 1.0):
        rows = [r for r in diagnostic if abs(float(r["theta"]) - theta) < 1e-12]
        t = np.array([float(r["elapsed_min"]) for r in rows])
        style = "-" if theta == 1 else "--"
        axes[0].plot(t, [float(r["A_mol_L"]) for r in rows], style,
                     label=f"θ={theta:g}, discrete A")
    rows = [r for r in diagnostic if abs(float(r["theta"]) - .5) < 1e-12]
    axes[0].plot([float(r["elapsed_min"]) for r in rows],
                 [float(r["exact_A_mol_L"]) for r in rows], color="black", lw=2,
                 label="Exact A")
    axes[0].axhline(0, color="black", lw=.7)
    axes[0].set_xlabel("Elapsed time (min)")
    axes[0].set_ylabel("A concentration (mol/L)")
    axes[0].set_title("Unresolved fast reaction: h = 0.01 min, k_f h = 10")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=.2)

    species_colors = {"A": "#2878b5", "B": "#ed8b00", "C": "#2b9348"}
    for schedule, style in (("early_sampled", "-"), ("late_sampled", "--")):
        rows = [r for r in trajectories if r["schedule"] == schedule]
        t = 60.0 * np.array([float(r["elapsed_min"]) for r in rows])
        for species_name in ("A", "B", "C"):
            axes[1].plot(t, [float(r[f"theta_{species_name}"]) for r in rows], style,
                         color=species_colors[species_name], alpha=.85,
                         label=f"{species_name}, {schedule.replace('_', ' ')}")
            if schedule == "early_sampled":
                axes[1].plot(t, [float(r[f"exact_{species_name}"]) for r in rows], ":",
                             color=species_colors[species_name], lw=1.3,
                             label=f"{species_name}, exact")
    product_data = [r for r in observations if r["schedule"] == "early_sampled"]
    axes[1].scatter([60.0 * float(r["elapsed_min"]) for r in product_data],
                    [float(r["observed_C_mol_L"]) for r in product_data],
                    color="#6930c3", marker="x", s=30, label="Product measurements", zorder=4)

    qualification_colors = {"early_sampled": "#6a3d9a", "late_sampled": "#b23a48"}
    qualification_styles = {"early_sampled": "--", "late_sampled": "-."}
    for schedule in ("early_sampled", "late_sampled"):
        by_goal = {r["goal"]: r for r in qualification if r["schedule"] == schedule}
        for first, second in (("A", "B"),):
            row_first, row_second = by_goal[first], by_goal[second]
            if row_first["decision"] != "qualified" or row_second["decision"] != "qualified":
                continue
            time_first = 60.0 * float(row_first["persistent_elapsed_min"])
            time_second = 60.0 * float(row_second["persistent_elapsed_min"])
            if abs(time_first - time_second) < 1e-10:
                axes[1].axvline(
                    time_first, color=qualification_colors[schedule],
                    linestyle=qualification_styles[schedule], linewidth=1.6,
                    label=f"A,B {schedule.replace('_', ' ')} qualification: {time_first:.4f} s",
                    zorder=3)
            else:
                for goal, row in ((first, row_first), (second, row_second)):
                    time_s = 60.0 * float(row["persistent_elapsed_min"])
                    axes[1].axvline(
                        time_s, color=qualification_colors[schedule],
                        linestyle=qualification_styles[schedule], linewidth=1.4,
                        label=f"{goal} {schedule.replace('_', ' ')} qualification: {time_s:.4f} s",
                        zorder=3)
    axes[1].set_xlim(0, 1.5)
    axes[1].set_xlabel("Elapsed time (s)")
    axes[1].set_ylabel("Concentration (mol/L)")
    axes[1].set_title("State reconstruction; C qualifies from $t_0$")
    axes[1].grid(alpha=.2)
    fig.suptitle("Stiff consecutive reaction A → B → C")
    fig.legend(*axes[1].get_legend_handles_labels(), loc="lower center",
               bbox_to_anchor=(.5, .005), ncol=4, fontsize=7.5, frameon=False)
    fig.tight_layout(rect=(0, .19, 1, .94))
    fig.savefig(root / "stiffness_and_reconstruction.png", dpi=180)
    plt.close(fig)

    error_fields = ["state_error", "tangent_error", "goal_adjoint_error", "initial_fit_error",
                    "objective_error", "gradient_error", "hessian_error", "fitted_final_state_error"]
    labels = ["State", "Tangent", "Goal adjoint", "Initial fit", "Objective at probe",
              "Gradient", "Hessian", "Fitted final state"]
    fig, axes = plt.subplots(2, 4, figsize=(15, 7))
    colors = {0.5: "#2455a4", 0.65: "#c44536"}
    for ax, field, label in zip(axes.flat, error_fields, labels):
        for theta in (0.5, 0.65):
            rows = [r for r in convergence if abs(float(r["theta"]) - theta) < 1e-12]
            ax.loglog([float(r["max_step_min"]) for r in rows],
                      [float(r[field]) for r in rows], "o-", ms=3,
                      color=colors[theta], label=f"θ={theta:g}")
        ax.set_title(label)
        ax.set_xlabel("Maximum time step (min)")
        ax.grid(True, which="both", alpha=.2)
    axes[0, 0].set_ylabel("Error")
    axes[1, 0].set_ylabel("Error")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Temporal convergence of the discrete state-finding system")
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(root / "convergence.png", dpi=180)
    plt.close(fig)

    quantity_labels = ["state", "tangent", "goal adjoint", "initial fit", "objective",
                       "gradient", "Hessian", "final state"]
    by_theta = {}
    for theta in (.5, .65):
        rows = [r for r in orders if abs(float(r["theta"]) - theta) < 1e-12]
        by_theta[theta] = rows[-8:]
    matrix = np.array([[float(r["observed_order"]) for r in by_theta[theta]]
                       for theta in (.5, .65)])
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.6), sharey=True)
    for i, (ax, theta, expected) in enumerate(zip(axes, (.5, .65), (2.0, 1.0))):
        ax.bar(np.arange(8), matrix[i], color=colors[theta])
        ax.axhline(expected, color="black", ls="--", label=f"Expected {expected:g}")
        ax.set_xticks(np.arange(8), quantity_labels, rotation=40, ha="right")
        ax.set_ylim(0, max(2.3, expected + .3))
        ax.set_title(f"θ={theta:g}")
        ax.set_ylabel("Last-pair observed order")
        ax.legend(frameon=False)
        for j, value in enumerate(matrix[i]):
            ax.text(j, value + .035, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    fig.suptitle("Measured convergence order across the workflow")
    fig.tight_layout()
    fig.savefig(root / "convergence_orders.png", dpi=180)
    plt.close(fig)

    obsolete_plot = root / "qualification_times.png"
    if obsolete_plot.exists():
        obsolete_plot.unlink()

    print(f"Plotted C++ stiff-reaction outputs in {root}; Python performs no model or qualification calculations.")


if __name__ == "__main__":
    main()
