#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Render retained C++ outputs for the winding/housing coupling study."""
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
    parser.add_argument("--results", default="results/weak_coupling", type=Path)
    args = parser.parse_args()
    root = args.results
    root.mkdir(parents=True, exist_ok=True)
    cases = read_csv(root / "summary.csv")
    trajectories = read_csv(root / "trajectories.csv")
    observations = read_csv(root / "observations.csv")
    envelopes = read_csv(root / "qualification_envelopes.csv")
    decisions = read_csv(root / "qualification_decisions.csv")
    convergence = read_csv(root / "convergence.csv")
    costates = read_csv(root / "costates.csv")

    colors = {"weak": "#b23a48", "moderate": "#e09f3e", "strong": "#227c9d"}
    goals = ["winding_temperature", "housing_temperature", "winding_housing_difference", "stored_heat"]
    goal_labels = ["Winding temperature", "Housing temperature", "Winding − housing", "Stored heat deviation"]

    fig, axes = plt.subplots(3, 2, figsize=(13, 11), sharex=True)
    for row, case in enumerate(cases):
        name = case["case"]
        path = [r for r in trajectories if r["case"] == name]
        obs = [r for r in observations if r["case"] == name]
        t = np.array([float(r["elapsed_min"]) for r in path])
        axw, axh = axes[row]
        axw.plot(t, [float(r["truth_exact_winding_C"]) for r in path], color="#222", lw=2.2, label="Exact truth")
        axw.plot(t, [float(r["optimized_theta_winding_C"]) for r in path], color=colors[name], ls="--", lw=1.7, label="Discrete reconstruction")
        axw.fill_between(t,
                         [float(r["optimized_theta_winding_C"]) - float(r["winding_numerical_allowance_K"]) for r in path],
                         [float(r["optimized_theta_winding_C"]) + float(r["winding_numerical_allowance_K"]) for r in path],
                         color=colors[name], alpha=.14, label="Numerical allowance")
        axh.plot(t, [float(r["truth_exact_housing_C"]) for r in path], color="#222", lw=2.2, label="Exact truth")
        axh.plot(t, [float(r["optimized_theta_housing_C"]) for r in path], color=colors[name], ls="--", lw=1.7, label="Discrete reconstruction")
        axh.fill_between(t,
                         [float(r["optimized_theta_housing_C"]) - float(r["housing_numerical_allowance_K"]) for r in path],
                         [float(r["optimized_theta_housing_C"]) + float(r["housing_numerical_allowance_K"]) for r in path],
                         color=colors[name], alpha=.14, label="Numerical allowance")
        axh.errorbar([float(r["elapsed_min"]) for r in obs], [float(r["observed_housing_C"]) for r in obs],
                     yerr=[float(r["joint_noise_box_halfwidth_K"]) for r in obs], fmt="o", ms=4,
                     color="#6930c3", capsize=2, label="Housing observations; joint 95% box")
        axw.set_ylabel("Winding temperature (°C)")
        axh.set_ylabel("Housing temperature (°C)")
        axw.set_title(f"{name.capitalize()} contact, K = {float(case['contact_kW_per_K']):g} kW/K")
        if row == 2:
            axw.set_xlabel("Elapsed time (min)")
            axh.set_xlabel("Elapsed time (min)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Hidden winding temperature reconstructed from housing measurements", fontsize=15)
    fig.tight_layout(rect=(0, .06, 1, .96))
    fig.savefig(root / "reconstruction.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for ax, case in zip(axes, cases):
        name = case["case"]
        for goal, label in zip(goals, goal_labels):
            rows = [r for r in envelopes if r["case"] == name and r["goal"] == goal]
            t = np.array([float(r["elapsed_min"]) for r in rows])
            ratio = np.array([float(r["bound"]) / float(r["tolerance"]) for r in rows])
            ax.plot(t, ratio, lw=1.6, label=label)
        ax.axhline(1.0, color="black", lw=1, ls=":")
        for row in decisions:
            if row["case"] == name and row["persistent_elapsed_min"] != "nan":
                ax.axvline(float(row["persistent_elapsed_min"]), color=f"C{goals.index(row['goal'])}",
                           lw=.8, alpha=.3)
        ax.set_yscale("log")
        ax.set_ylabel("Bound / tolerance")
        ax.set_title(f"{name.capitalize()} contact")
        ax.grid(True, which="both", alpha=.2)
    axes[-1].set_xlabel("Elapsed time (min); line marks each goal's persistent start")
    axes[0].legend(ncol=2, frameon=False)
    fig.suptitle("Whole-interval uncertainty bounds and persistent qualification", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(root / "qualification.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(4, 4, figsize=(15, 12))
    error_fields = ["state_error_inf_K", "tangent_error_inf", "goal_adjoint_error_inf",
                    "initial_fit_error_K", "objective_error_at_truth", "gradient_error",
                    "hessian_error", "fitted_final_state_error_K"]
    titles = ["State", "Tangent", "Goal adjoint", "Recovered initial state", "Objective", "Gradient", "Hessian", "Final reconstructed state"]
    for theta_index, theta in enumerate((0.5, 0.65)):
        for index, (field, title) in enumerate(zip(error_fields, titles)):
            row = 2 * theta_index + index // 4
            col = index % 4
            ax = axes[row, col]
            for case in cases:
                name = case["case"]
                items = [r for r in convergence if r["case"] == name and abs(float(r["theta"]) - theta) < 1e-12]
                h = np.array([float(r["max_step_min"]) for r in items])
                error = np.array([float(r[field]) for r in items])
                ax.loglog(h, error, "o-", ms=3, lw=1.1, color=colors[name], label=name)
            ax.set_title(title)
            ax.set_xlabel("Maximum step (min)")
            if col == 0:
                ax.set_ylabel(f"θ = {theta:g}\nabsolute error")
            ax.grid(True, which="both", alpha=.18)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Observed convergence across the state-finding and discrete optimality workflow", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(root / "convergence.png", dpi=180)
    plt.close(fig)

    order_rows = read_csv(root / "orders.csv")
    quantities = ["state", "tangent", "goal_adjoint", "initial_fit", "objective",
                  "gradient", "hessian", "fitted_final_state"]
    quantity_labels = ["State", "Tangent", "Goal adjoint", "Initial fit", "Objective",
                       "Gradient", "Hessian", "Final state"]
    case_names = ["weak", "moderate", "strong"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8), sharey=True)
    for ax, theta, expected in zip(axes, (0.5, 0.65), (2.0, 1.0)):
        matrix = np.empty((len(case_names), len(quantities)))
        for i, case_name in enumerate(case_names):
            for j, quantity in enumerate(quantities):
                row = next(r for r in order_rows
                           if abs(float(r["contact_kW_per_K"]) - float(next(
                               c["contact_kW_per_K"] for c in cases if c["case"] == case_name))) < 1e-12
                           and abs(float(r["theta"]) - theta) < 1e-12
                           and r["quantity"] == quantity)
                matrix[i, j] = float(row["observed_order"])
        deviation = matrix - expected
        image = ax.imshow(deviation, cmap="coolwarm", vmin=-0.22, vmax=0.22, aspect="auto")
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                ax.text(j, i, f"{matrix[i, j]:.3f}", ha="center", va="center", fontsize=8)
        ax.set_xticks(range(len(quantity_labels)), quantity_labels, rotation=35, ha="right")
        ax.set_yticks(range(len(case_names)), [s.capitalize() for s in case_names])
        ax.set_title(f"θ = {theta:g}: expected order {expected:g}")
        ax.set_xlabel("Quantity")
        ax.set_xticks(np.arange(-.5, len(quantity_labels), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(case_names), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.2)
        ax.tick_params(which="minor", bottom=False, left=False)
    axes[0].set_ylabel("Thermal contact case")
    colorbar = fig.colorbar(image, ax=axes, shrink=.82, pad=.03)
    colorbar.set_label("Observed order − expected order")
    fig.suptitle("Measured convergence rates from the final refinement pair")
    fig.savefig(root / "convergence_orders.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    contact = np.array([float(r["contact_kW_per_K"]) for r in cases])
    half = np.array([float(r["slow_half_life_min"]) for r in cases])
    cond = np.array([float(r["scaled_information_condition"]) for r in cases])
    axes[0].semilogx(contact, half, "o-", color="#6930c3", lw=2)
    axes[0].set_xlabel("Winding-to-housing conductance K (kW/K)")
    axes[0].set_ylabel("Slow-mode half-life (min)")
    axes[0].grid(True, which="both", alpha=.25)
    axes[1].loglog(contact, cond, "s-", color="#d1495b", lw=2)
    axes[1].set_xlabel("Winding-to-housing conductance K (kW/K)")
    axes[1].set_ylabel("Scaled information-matrix condition number")
    axes[1].grid(True, which="both", alpha=.25)
    fig.suptitle("Thermal timescale separation and information about the hidden winding state")
    fig.tight_layout()
    fig.savefig(root / "modes_and_information.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for ax, case in zip(axes, cases):
        rows = [r for r in costates if r["case"] == case["case"]]
        t = [float(r["elapsed_min"]) for r in rows]
        ax.plot(t, [float(r["lambda_winding"]) for r in rows], label="Winding costate")
        ax.plot(t, [float(r["lambda_housing"]) for r in rows], label="Housing costate")
        ax.set_ylabel("Costate (K⁻¹)")
        ax.set_title(f"{case['case'].capitalize()} contact")
        ax.grid(True, alpha=.2)
    axes[-1].set_xlabel("Elapsed time (min)")
    axes[0].legend(frameon=False)
    fig.suptitle("Discrete costate: how each state contributes to future housing-data mismatch")
    fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(root / "costates.png", dpi=180)
    plt.close(fig)
    print(f"Plotted C++ weak-coupling outputs in {root}, including convergence_orders.png; no model or qualification logic runs in Python.")


if __name__ == "__main__":
    main()
