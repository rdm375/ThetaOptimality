#!/usr/bin/env python3
"""Plot C++ predator-prey reconstruction and convergence CSV outputs only."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "predator_prey"


def read_csv(name):
    with (RESULTS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    reconstruction = read_csv("reconstruction.csv")
    observations = read_csv("observations.csv")
    t = [float(row["time"]) for row in reconstruction]
    cases = [
        ("initial only", "prey_initial_only", "predator_initial_only", "tab:blue"),
        ("rates only", "prey_parameters_only", "predator_parameters_only", "tab:orange"),
        ("joint", "prey_joint", "predator_joint", "tab:green"),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    for axis, species, sd_key, title in [
        (axes[0], "prey", "prey_sd", "Prey population"),
        (axes[1], "predator", "predator_sd", "Predator population"),
    ]:
        axis.plot(t, [float(r[f"{species}_reference"]) for r in reconstruction],
                  "k--", linewidth=2.2, label="continuous reference")
        for label, prey_key, predator_key, color in cases:
            key = prey_key if species == "prey" else predator_key
            axis.plot(t, [float(r[key]) for r in reconstruction], color=color,
                      label=f"theta reconstruction: {label}")
        obs_key = f"{species}_observed"
        axis.errorbar([float(r["time"]) for r in observations],
                      [float(r[obs_key]) for r in observations],
                      yerr=[float(r[sd_key]) for r in observations],
                      fmt="ko", markersize=3.5, capsize=2, label="synthetic assays")
        axis.set_ylabel("normalized population")
        axis.set_title(title)
        axis.grid(True, alpha=.25)
    axes[1].set_xlabel("time")
    axes[0].legend(ncol=2, fontsize=8)
    fig.suptitle("Predator–prey reconstruction: known state, known rates, and joint recovery")
    fig.savefig(RESULTS / "state_reconstruction.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8, 7), constrained_layout=True)
    axis.plot([float(r["prey_reference"]) for r in reconstruction],
              [float(r["predator_reference"]) for r in reconstruction],
              "k--", linewidth=2.2, label="continuous reference")
    for label, prey_key, predator_key, color in cases:
        axis.plot([float(r[prey_key]) for r in reconstruction],
                  [float(r[predator_key]) for r in reconstruction],
                  color=color, label=f"theta reconstruction: {label}")
    axis.scatter([float(r["prey_observed"]) for r in observations],
                 [float(r["predator_observed"]) for r in observations],
                 s=20, color="black", marker="x", label="assay pairs")
    axis.set_xlabel("prey population")
    axis.set_ylabel("predator population")
    axis.set_title("Phase portrait and the effect of which quantities are estimated")
    axis.grid(True, alpha=.25)
    axis.legend()
    fig.savefig(RESULTS / "phase_portrait.png", dpi=180)
    plt.close(fig)

    rows = read_csv("convergence.csv")
    labels = {"initial_state": ("initial populations", "tab:blue"),
              "parameters": ("rates", "tab:orange"),
              "joint": ("joint", "tab:green")}
    styles = {"0.5": "-o", "0.65": "--s"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    panels = [("state_error", "continuous trajectory error"),
              ("fit_error", "change in estimated controls vs finest grid"),
              ("gradient_error", "adjoint gradient finite-difference error"),
              ("hessian_error", "incremental-adjoint Hessian finite-difference error")]
    for axis, (key, ylabel) in zip(axes.flat, panels):
        for case, (case_label, color) in labels.items():
            for theta, style in styles.items():
                selected = [r for r in rows if r["case"] == case and r["theta"] == theta]
                x = [float(r["h"]) for r in selected]
                y = [float(r[key]) for r in selected]
                if key == "fit_error":
                    pairs = [(a, b) for a, b in zip(x, y) if b > 0]
                    x, y = zip(*pairs)
                axis.loglog(x, y, style, color=color, label=f"{case_label}, theta={theta}")
        axis.set_xlabel("step size h")
        axis.set_ylabel(ylabel)
        axis.grid(True, which="both", alpha=.25)
        axis.legend(fontsize=7)
    fig.suptitle("Theta refinement and derivative checks for all three inverse problems")
    fig.savefig(RESULTS / "convergence_and_derivatives.png", dpi=180)
    plt.close(fig)

    print("Plotted retained C++ predator-prey outputs; Python performs no simulation, fitting, or order estimation.")


if __name__ == "__main__":
    main()
