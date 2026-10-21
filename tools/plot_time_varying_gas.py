#!/usr/bin/env python3
"""Plot C++ outputs only; no simulation, fitting, or qualification runs here."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "time_varying_gas"


def read(name):
    with (OUT / name).open(newline="") as f:
        return list(csv.DictReader(f))


fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
rows = read("state_theta_050.csv")
t = [float(r["time"]) for r in rows]
for ax, key, title in zip(
    axes,
    ("ethane", "nitrogen", "methane"),
    ("Ethane mole fraction", "Nitrogen mole fraction", "Methane mole fraction"),
):
    # Draw the reconstruction first so the analytic truth remains visible
    # when both curves nearly coincide.
    ax.plot(t, [float(r[f"{key}_fit"]) for r in rows], color="tab:blue",
            lw=2.0, alpha=.8, label="Reconstruction", zorder=2)
    ax.plot(t, [float(r[f"{key}_true"]) for r in rows], "k:",
            lw=2.2, marker="o", ms=2.2, markevery=20, label="Analytical truth",
            zorder=3)
    ax.set_ylabel("Mole fraction")
    ax.set_title(title)
    ax.grid(alpha=.25)
axes[-1].set_xlabel("Elapsed time (min)")
axes[0].legend()
fig.suptitle("Well-mixed vessel composition, midpoint reconstruction")
fig.tight_layout()
fig.savefig(OUT / "composition_reconstruction.png", dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 4.5))
for filename, label, style in (
    ("state_theta_050.csv", r"$\theta=0.5$", "-"),
    ("state_theta_065.csv", r"$\theta=0.65$", "--"),
):
    rows = read(filename)
    t = [float(r["time"]) for r in rows]
    ax.plot(t, [float(r["feed_ethane_true"]) for r in rows], "k:", lw=1.7)
    ax.plot(t, [float(r["feed_ethane_fit"]) for r in rows], style, label=f"Fitted feed {label}")
ax.plot([], [], "k:", label="True feed schedule")
ax.set(xlabel="Elapsed time (min)", ylabel="Ethane mole fraction in feed",
       title="Recovered time-varying feed composition")
ax.grid(alpha=.25)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "feed_reconstruction.png", dpi=180)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
rows = read("convergence.csv")
for theta, style in ((.5, "o-"), (.65, "s--")):
    subset = [r for r in rows if abs(float(r["theta"]) - theta) < 1e-9]
    axes[0].loglog([int(r["N"]) for r in subset],
                   [float(r["state_error"]) for r in subset], style,
                   label=fr"$\theta={theta}$")
    axes[1].loglog([int(r["N"]) for r in subset],
                   [float(r["parameter_error"]) for r in subset], style,
                   label=fr"$\theta={theta}$")
axes[0].set(title="Analytical state error", xlabel="Time steps", ylabel="Max state error")
axes[1].set(title="Joint initial/feed coefficient error", xlabel="Time steps", ylabel="Euclidean error")
for ax in axes:
    ax.grid(True, which="both", alpha=.25)
    ax.legend()
fig.suptitle("Time refinement for the dynamic gas reconstruction")
fig.tight_layout()
fig.savefig(OUT / "convergence_and_parameter_error.png", dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 4.2))
rows = read("density_residuals.csv")
for theta, marker in ((.5, "o"), (.65, "s")):
    subset = [r for r in rows if abs(float(r["theta"]) - theta) < 1e-9]
    ax.plot([float(r["time"]) for r in subset],
            [float(r["residual"]) / float(r["sigma"]) for r in subset],
            marker + "-", ms=3, label=fr"$\theta={theta}$")
ax.axhline(0, color="black", lw=.8)
ax.set(xlabel="Elapsed time (min)", ylabel=r"Density residual / $\sigma$",
       title="Peng–Robinson observation residuals")
ax.grid(alpha=.25)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "density_residuals.png", dpi=180)
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 4.5))
rows = read("information_spectrum.csv")
for theta, marker in ((.5, "o-"), (.65, "s--")):
    row = next(r for r in rows if abs(float(r["theta"]) - theta) < 1e-9)
    eigenvalues = [float(row[f"eigenvalue_{j}"]) for j in range(5)]
    ax.semilogy(range(1, 6), eigenvalues, marker, label=fr"$\theta={theta}$")
ax.set(xlabel="Information mode (largest to smallest)", ylabel="Eigenvalue",
       title="Local information in initial composition and feed knots")
ax.grid(True, which="both", alpha=.25)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "information_spectrum.png", dpi=180)
plt.close(fig)
print("Plotted C++ dynamic gas outputs; Python performs no thermodynamic evaluation, fitting, or convergence analysis.")
