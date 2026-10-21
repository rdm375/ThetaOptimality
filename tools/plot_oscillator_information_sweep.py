#!/usr/bin/env python3
"""Render retained C++ weak-information sweep data; performs no inference."""
from pathlib import Path
import csv
import math
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "oscillator_information_sweep"

def read_csv(name):
    with (OUT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))

rows = read_csv("sweep.csv")
alpha = [float(r["alpha"]) for r in rows]
small = [max(float(r["analytic_lambda_small"]), 1e-25) for r in rows]
discrete = [max(float(r["discrete_lambda_small"]), 1e-25) for r in rows]
vel_width = [float(r["v0_bounded_halfwidth_m_s"]) for r in rows]
vel_sd = [float(r["v0_local_sd_m_s"]) for r in rows]
vfit = [float(r["v0_fit_m_s"]) for r in rows]
grad = [max(float(r["gradient_norm"]), 1e-18) for r in rows]

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
ax = axes[0, 0]
ax.semilogy(alpha, small, "o-", label="exact physical information")
ax.semilogy(alpha, discrete, "s--", label="midpoint discrete Hessian")
ax.set(xlabel="design deformation $\\alpha$", ylabel="small scaled information eigenvalue",
       title="Information vanishes at the blind schedule")
ax.grid(True, which="both", alpha=.25); ax.legend()
ax.annotate("exact rank 1", (1, small[-1]), xytext=(.76, 8e-23),
            arrowprops={"arrowstyle": "->"})

ax = axes[0, 1]
ax.semilogy(alpha, vel_width, "o-", label="bounded-data half-width")
ax.semilogy(alpha, vel_sd, "s--", label="local Gaussian standard deviation")
ax.axhline(.002, color="black", linestyle=":", label="velocity tolerance")
ax.set(xlabel="design deformation $\\alpha$", ylabel="initial-velocity uncertainty (m/s)",
       title="Practical uncertainty grows as timing degrades")
ax.grid(True, which="both", alpha=.25); ax.legend()

ax = axes[1, 0]
ax.plot(alpha, vfit, "o-", label="discrete fitted initial velocity")
ax.axhline(.12, color="black", linestyle="--", label="synthetic truth")
ax.set(xlabel="design deformation $\\alpha$", ylabel="initial velocity (m/s)",
       title="A stationary fit can become physically arbitrary")
ax.grid(True, alpha=.25); ax.legend()
ax2 = ax.twinx(); ax2.semilogy(alpha, grad, color="tab:red", linestyle=":", label="optimizer gradient norm")
ax2.set_ylabel("discrete objective gradient norm", color="tab:red")
ax2.tick_params(axis="y", labelcolor="tab:red")
ax.annotate("large phase-leakage fit", (1, vfit[-1]), xytext=(.67, -1.7),
            arrowprops={"arrowstyle": "->"})

ax = axes[1, 1]
for key, label, color in (("q0_exact_qualified_t0", "displacement", "tab:blue"),
                          ("v0_exact_qualified_t0", "velocity", "tab:orange"),
                          ("support_force_exact_qualified_t0", "support force", "tab:green")):
    values = [int(r[key]) for r in rows]
    ax.step(alpha, values, where="mid", label=label, color=color)
ax.set(xlabel="design deformation $\\alpha$", ylabel="exact-model qualified at $t_0$",
       yticks=[0, 1], yticklabels=["no", "yes"], title="Bounded-data goal qualification")
ax.grid(True, alpha=.25); ax.legend(loc="best")
fig.suptitle("Oscillator weak-information timing sweep\nExact-model qualification excludes discretization error", fontsize=14)
fig.savefig(OUT / "weak_information_sweep.png", dpi=180)
plt.close(fig)

schedule_rows = read_csv("sampling_schedules.csv")
fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
for a in (0., .25, .5, .75, 1.):
    selected = [r for r in schedule_rows if math.isclose(float(r["alpha"]), a, abs_tol=1e-12)]
    times = [float(r["elapsed_s"]) for r in selected]
    ax.scatter(times, [a] * len(times), s=55, label=f"$\\alpha={a:g}$")
ax.set(xlabel="elapsed sample time (s)", ylabel="design deformation $\\alpha$",
       title="Four displacement samples slide toward half-periods", yticks=[0, .25, .5, .75, 1])
ax.grid(True, alpha=.25); ax.legend(ncol=5, loc="upper center")
fig.savefig(OUT / "sampling_schedule_sweep.png", dpi=180)
plt.close(fig)
print("Plotted the C++ weak-information sweep; Python performs no fitting or qualification.")
