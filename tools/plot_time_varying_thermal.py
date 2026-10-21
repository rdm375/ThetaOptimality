#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot retained C++ results; this script does no simulation or inference."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "time_varying_thermal"

def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))

cases = (("050", r"$\theta=0.5$", "tab:blue", "-"),
         ("065", r"$\theta=0.65$", "tab:orange", "--"))

# State trajectories and the reconstructed, time-varying heater.
fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
for suffix, label, color, style in cases:
    rows = read_csv(RESULTS / f"state_theta_{suffix}.csv")
    t = [float(row["time"]) for row in rows]
    for key, name, state_color in (("fluid", "Fluid", "tab:red"), ("wall", "Wall", "tab:blue")):
        axes[0].plot(t, [float(row[f"{key}_fit"]) for row in rows], style, color=state_color,
                     label=f"{name}, fitted {label}", alpha=.9)
    axes[1].plot(t, [float(row["heater_fit"]) for row in rows], style, color=color,
                 label=f"Fitted heater {label}")
    if suffix == "050":
        axes[0].plot(t, [float(row["fluid_reference"]) for row in rows], ":", color="tab:red", lw=1.5,
                     label="Analytical fluid state")
        axes[0].plot(t, [float(row["wall_reference"]) for row in rows], ":", color="tab:blue", lw=1.5,
                     label="Analytical wall state")
        axes[1].plot(t, [float(row["heater_truth"]) for row in rows], ":", color="black", lw=1.5,
                     label="True heater schedule")
observed = read_csv(RESULTS / "observation_residuals.csv")
observed_midpoint = [r for r in observed if float(r["theta"]) == .5]
axes[0].scatter([float(r["time"]) for r in observed_midpoint], [float(r["observation"]) for r in observed_midpoint],
                s=16, color="black", zorder=5, label="Fluid observations")
axes[0].set_ylabel("Temperature (°C)")
axes[1].set_ylabel("Heater input (scaled)")
axes[1].set_xlabel("Elapsed time (min)")
axes[0].legend(ncol=2, fontsize=7)
axes[1].legend(ncol=2, fontsize=8)
axes[0].grid(alpha=.25); axes[1].grid(alpha=.25)
fig.suptitle("State reconstruction and heater control recovery")
fig.savefig(RESULTS / "reconstruction.png", dpi=180)

# Mesh-refinement errors. The joint coefficient error shows that a close state
# trajectory does not imply equally accurate heater-knot estimates.
refinement = read_csv(RESULTS / "convergence.csv")
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
for suffix, label, color, _ in cases:
    theta = .5 if suffix == "050" else .65
    rows = [r for r in refinement if float(r["theta"]) == theta]
    h = [10 / int(r["steps"]) for r in rows]
    axes[0].loglog(h, [float(r["trajectory_error"]) for r in rows], "o-", color=color,
                   label=f"State, {label}")
    axes[1].loglog(h, [float(r["fit_error"]) for r in rows], "o-", color=color,
                   label=f"All fitted coefficients, {label}")
for ax in axes:
    ax.invert_xaxis(); ax.set_xlabel("Step size $h$ (min)"); ax.grid(True, which="both", alpha=.25)
axes[0].set_ylabel("Maximum state error (°C)")
axes[1].set_ylabel("Maximum coefficient error")
axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
fig.suptitle("State and joint-estimate convergence under mesh refinement")
fig.savefig(RESULTS / "convergence.png", dpi=180)

# Each fitted variable is shown separately against its exact synthetic value.
parameters = read_csv(RESULTS / "parameter_refinement.csv")
names = ["fluid_initial", "wall_initial", "heater_knot_0", "heater_knot_5", "heater_knot_10"]
labels = ["Initial fluid temperature (°C)", "Initial wall temperature (°C)",
          "Heater at 0 min (scaled)", "Heater at 5 min (scaled)", "Heater at 10 min (scaled)"]
fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
for j, (name, ylabel) in enumerate(zip(names, labels)):
    ax = axes.flat[j]
    subset = [r for r in parameters if r["parameter"] == name]
    truth = float(subset[0]["truth"])
    ax.axhline(truth, color="black", linestyle=":", label="Analytical truth")
    for suffix, label, color, _ in cases:
        theta = .5 if suffix == "050" else .65
        rows = [r for r in subset if float(r["theta"]) == theta]
        ax.plot([int(r["steps"]) for r in rows], [float(r["estimate"]) for r in rows],
                "o-", color=color, label=label)
    ax.set_xscale("log", base=2); ax.set_xlabel("Time steps $N$"); ax.set_ylabel(ylabel)
    ax.grid(True, alpha=.25)
    if j == 0: ax.legend(fontsize=8)
axes.flat[-1].axis("off")
fig.suptitle("Recovered initial temperatures and heater knots")
fig.savefig(RESULTS / "parameter_recovery.png", dpi=180)

# Residuals relative to the stated measurement scale.
fig, ax = plt.subplots(figsize=(8.5, 4.2), constrained_layout=True)
for suffix, label, color, _ in cases:
    theta = .5 if suffix == "050" else .65
    rows = [r for r in observed if float(r["theta"]) == theta]
    ax.plot([float(r["time"]) for r in rows],
            [float(r["residual"]) / float(r["measurement_scale"]) for r in rows],
            "o-", color=color, label=label)
ax.axhline(0., color="black", lw=.8)
ax.set_xlabel("Elapsed time (min)"); ax.set_ylabel("Residual / measurement scale")
ax.grid(alpha=.25); ax.legend()
ax.set_title("One-thermometer observation residuals at $N=640$")
fig.savefig(RESULTS / "observation_residuals.png", dpi=180)

# Local information eigenvalues diagnose how well the five coefficients are
# separated by the sampled fluid-temperature record.
spectrum = read_csv(RESULTS / "information_spectrum.csv")
fig, ax = plt.subplots(figsize=(8.5, 4.2), constrained_layout=True)
indices = list(range(5)); width = .36
for k, (suffix, label, color, _) in enumerate(cases):
    theta = .5 if suffix == "050" else .65
    row = next(r for r in spectrum if float(r["theta"]) == theta)
    vals = [float(row[f"eigenvalue_{j}"]) for j in indices]
    ax.bar([j + (k-.5)*width for j in indices], vals, width, color=color, label=label)
ax.set_yscale("log"); ax.set_xticks(indices, [f"$\\lambda_{j+1}$" for j in indices])
ax.set_ylabel("Eigenvalue of scaled observation information")
ax.set_xlabel("Information mode (largest to smallest)")
ax.grid(True, axis="y", which="both", alpha=.25); ax.legend()
ax.set_title("Local information spectrum for the joint state-control fit")
fig.savefig(RESULTS / "information_spectrum.png", dpi=180)

# Local linearized uncertainty under the stated independent Gaussian sensor
# scale. These are diagnostics, not the physical engineering qualification.
uncertainty = read_csv(RESULTS / "local_parameter_uncertainty.csv")
fig, axes = plt.subplots(1, 5, figsize=(14, 3.8), constrained_layout=True)
for j, (name, ylabel) in enumerate(zip(names, labels)):
    ax = axes[j]
    subset = [r for r in uncertainty if r["parameter"] == name]
    vals = [float(next(r for r in subset if float(r["theta"]) == theta)["local_standard_deviation"])
            for theta in (.5, .65)]
    ax.bar([0, 1], vals, color=["tab:blue", "tab:orange"])
    ax.set_xticks([0, 1], [r"$\theta=.5$", r"$\theta=.65$"])
    ax.set_title(ylabel, fontsize=8); ax.set_ylabel("Local SD")
    ax.grid(True, axis="y", alpha=.25)
fig.suptitle("Local parameter uncertainty from the one-sensor information matrix")
fig.savefig(RESULTS / "local_uncertainty.png", dpi=180)
print("Plotted state/control reconstruction, refinement, coefficient recovery, observation residuals, information spectrum, and local uncertainty.")
