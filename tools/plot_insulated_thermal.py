#!/usr/bin/env python3
"""Render retained C++ results; no forward, inverse, or qualification logic."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT = Path(__file__).resolve().parents[1] / 'results' / 'insulated_thermal'

def read(name):
    with (ROOT / name).open() as f:
        return list(csv.DictReader(f))

def numbers(rows, key):
    return np.array([float(r[key]) if r[key] else np.nan for r in rows])

def save(fig, name):
    fig.savefig(ROOT / (name + '.png'), dpi=170, bbox_inches='tight')
    plt.close(fig)

plt.rcParams.update({'axes.grid': True, 'grid.alpha': .2, 'font.size': 10})
paths = read('trajectories.csv')
cases = [('difference_only', 'Difference data | seed (24,24)', '#0072B2'),
         ('shifted_seed', 'Same data | seed (28,28)', '#D55E00'),
         ('added_absolute', 'Added absolute-temperature datum', '#009E73')]
first = [r for r in paths if r['case'] == 'difference_only']
t = numbers(first, 'elapsed_min')
fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
for ax, quantity, truthkey, unit in zip(axes.flat,
        ['t1', 't2', 'difference', 'mean'], ['true_t1', 'true_t2', None, 'true_mean'],
        ['Temperature 1 (°C)', 'Temperature 2 (°C)', 'Temperature difference (K)', 'Weighted mean temperature (°C)']):
    if truthkey:
        ax.plot(t, numbers(first, truthkey), color='.5', lw=3, label='Synthetic truth')
    for case, label, color in cases:
        rows = [r for r in paths if r['case'] == case]
        ax.plot(numbers(rows, 'elapsed_min'), numbers(rows, 'discrete_' + quantity), color=color, label=label)
        exactkey = 'exact_fit_' + quantity if quantity in ('t1', 't2') else 'exact_difference' if quantity == 'difference' else None
        if exactkey:
            ax.plot(numbers(rows, 'elapsed_min'), numbers(rows, exactkey), '--', color=color, lw=1)
    ax.set(xlabel='Elapsed time from t₀ (min)', ylabel=unit)
    if quantity == 'difference':
        obs = [r for r in read('observations.csv') if r['added'] == '0']
        ax.errorbar(numbers(obs, 'elapsed_min'), numbers(obs, 'value'), yerr=numbers(obs, 'scale'), fmt='ko', ms=4, label='Difference data ± bound')
axes[0, 0].legend(fontsize=8)
fig.suptitle('Insulated masses | discrete midpoint N=128; dashed curves are analytical fits\nDistinct absolute states fit the same difference record')
save(fig, 'reconstruction')

fig, ax = plt.subplots(figsize=(8, 7), layout='constrained')
ax.plot(numbers(first, 'true_t1'), numbers(first, 'true_t2'), color='.5', lw=3, label='Synthetic truth')
for case, label, color in cases:
    rows = [r for r in paths if r['case'] == case]
    ax.plot(numbers(rows, 'discrete_t1'), numbers(rows, 'discrete_t2'), color=color, label=label)
    ax.plot(numbers(rows, 'exact_fit_t1'), numbers(rows, 'exact_fit_t2'), '--', color=color)
    ax.scatter(float(rows[0]['discrete_t1']), float(rows[0]['discrete_t2']), color=color, marker='o')
    ax.scatter(float(rows[-1]['discrete_t1']), float(rows[-1]['discrete_t2']), color=color, marker='s')
ax.set(xlabel='Temperature 1 (°C)', ylabel='Temperature 2 (°C)', title='Phase paths | circle: initial; square: final\nDifferent conserved means remain separated')
ax.legend(fontsize=8)
save(fig, 'phase')

polygons = read('qualification_polygons.csv')
fig, axes = plt.subplots(1, 2, figsize=(12, 6), layout='constrained')
for ax, case, title in zip(axes, ['difference_refined', 'added_absolute'], ['Difference-only data: common-offset strip', 'Added absolute datum: local detail']):
    for which, color, style in [('outer', '#0072B2', '-'), ('analytical', '#D55E00', '--')]:
        rows = [r for r in polygons if r['case'] == case and r['set'] == which]
        rows = rows + rows[:1]
        ax.plot(numbers(rows, 't10'), numbers(rows, 't20'), style, color=color, label=which)
    ax.scatter([32], [20], c='k', marker='*', label='Truth', zorder=3)
    ax.set(xlabel='Initial temperature 1 (°C)', ylabel='Initial temperature 2 (°C)', title=title)
    ax.legend(fontsize=8)
save(fig, 'feasible_sets')

errors, orders = read('convergence.csv'), read('orders.csv')
quantities = [key for key in errors[0] if key not in ('added', 'theta', 'steps', 'h_min')]
fig, axes = plt.subplots(3, 4, figsize=(17, 12), layout='constrained')
for ax, key in zip(axes.flat, quantities):
    for added, theta, color, marker, style in [(0, .5, '#0072B2', 'o', '-'), (0, .65, '#D55E00', 's', '-'), (1, .5, '#009E73', '^', '--'), (1, .65, '#CC79A7', 'v', '--')]:
        rows = [r for r in errors if int(r['added']) == added and float(r['theta']) == theta]
        order = next(float(r['order']) for r in orders if int(r['added']) == added and float(r['theta']) == theta and r['quantity'] == key)
        ax.loglog(numbers(rows, 'h_min'), numbers(rows, key), marker=marker, linestyle=style, color=color, ms=4, label=f'{"Added" if added else "Difference"}, θ={theta}: {order:.3f}')
    ax.set(title=key.replace('_', ' '), xlabel='Step size (min)', ylabel='Error (stated scales)')
    ax.legend(fontsize=6)
fig.suptitle('Analytical convergence | representative-control error uses the same explicit gauge\nA converged representative does not establish an observable absolute state')
save(fig, 'convergence')

bounds, decisions = read('qualification_envelopes.csv'), read('qualification_decisions.csv')
goals = ['temperature1', 'temperature2', 'difference', 'weighted_mean']
titles = ['Temperature 1', 'Temperature 2', 'Difference', 'Weighted mean']
shown = [('difference_coarse', 'Difference data | θ=.5, N=32'), ('difference_refined', 'Difference data | θ=.5, N=512'), ('difference_theta65', 'Difference data | θ=.65, N=512'), ('added_absolute', 'Added absolute datum | θ=.5, N=512')]
fig, axes = plt.subplots(4, 4, figsize=(18, 14), layout='constrained')
for i, (case, label) in enumerate(shown):
    for j, goal in enumerate(goals):
        ax = axes[i, j]
        rows = [r for r in bounds if r['case'] == case and r['goal'] == goal]
        decision = next(r for r in decisions if r['case'] == case and r['goal'] == goal)
        ax.semilogy(numbers(rows, 'elapsed_min'), numbers(rows, 'whole_bound') / numbers(rows, 'tolerance'), color='#0072B2', label='Physical whole-interval bound')
        ax.semilogy(numbers(rows, 'elapsed_min'), numbers(rows, 'analytical_bound') / numbers(rows, 'tolerance'), '--', color='#D55E00', label='Analytical feasible-set reference')
        ax.axhline(1, color='k', lw=1)
        verdict = 'NOT QUALIFIED'
        color = '#B74425'
        if decision['qualified'] == '1':
            time = float(decision['persistent_min'])
            verdict = f'QUALIFIED from {time:.8g} min'
            color = '#14814B'
            ax.axvline(time, color=color, linestyle='--')
        ax.set_title(f'{titles[j]} | tolerance {float(decision["tolerance"]):g} {"K" if goal == "difference" else "°C"}\n{verdict}', color=color, fontsize=10)
        if j == 0:
            ax.set_ylabel(label + '\nError bound / tolerance', fontsize=9)
        if i == 3:
            ax.set_xlabel('Elapsed time (min)')
        if i == 0 and j == 0:
            ax.legend(fontsize=7)
fig.suptitle('Persistent engineering qualification | every remaining whole interval through 8 min\nBounded sensor errors and prior [10,40] °C × [10,40] °C; retrospective all-data decisions')
save(fig, 'qualification')

fig, ax = plt.subplots(figsize=(15, 5), layout='constrained')
ax.axis('off')
rows, colors = [], []
for case, label in shown:
    ds = [next(r for r in decisions if r['case'] == case and r['goal'] == goal) for goal in goals]
    row = [label] + [f'From {float(r["persistent_min"]):.8g} min' if r['qualified'] == '1' else 'NOT QUALIFIED' for r in ds]
    row.append(f'From {float(ds[0]["joint_min"]):.8g} min' if ds[0]['joint_qualified'] == '1' else 'NOT QUALIFIED')
    rows.append(row)
    colors.append(['#F0F2F3'] + ['#D7EBDF' if x.startswith('From') else '#F4DEDA' for x in row[1:]])
table = ax.table(cellText=rows, cellColours=colors, colLabels=['Experiment'] + titles + ['ALL GOALS'], cellLoc='center', loc='center', colWidths=[.28, .144, .144, .144, .144, .144])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 3)
ax.set_title('INSULATED THERMAL ENGINEERING QUALIFICATION\nDifference-only data cannot determine a persistent common temperature offset', fontsize=15, pad=20)
fig.text(.5, .02, 'Elapsed times; numerical error included. Successful validation includes correctly unqualified goals.', ha='center', fontsize=10)
save(fig, 'qualification_verdict')

rows = read('costates.csv')
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
for j, ax in enumerate(axes, 1):
    ax.plot(numbers(rows, 'elapsed_min'), numbers(rows, f'analytical{j}'), 'k--', lw=2, label='Analytical')
    ax.plot(numbers(rows, 'elapsed_min'), numbers(rows, f'discrete{j}'), color='#0072B2', label='Discrete cost-to-go')
    ax.set(xlabel='Elapsed time (min)', ylabel='Running-cost sensitivity (K min)', title=f'Running difference-cost costate, component {j}')
    ax.legend()
save(fig, 'costates')
print('Rendered C++ insulated-thermal outputs; Python performs no numerical inference or qualification.')
