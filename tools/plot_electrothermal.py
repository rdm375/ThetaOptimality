#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Render C++ results; no numerical model, optimization or qualification here."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
ROOT = Path(__file__).resolve().parents[1] / 'results' / 'electrothermal'

def read(name):
    with (ROOT / name).open() as f:
        return list(csv.DictReader(f))

def column(rows, key):
    return np.array([float(r[key]) if r[key] else np.nan for r in rows])

def save(fig, name):
    fig.savefig(ROOT / (name + '.png'), dpi=170, bbox_inches='tight')
    plt.close(fig)

plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': .2})
paths = read('trajectories.csv')
positive = [r for r in paths if r['branch'] == '1']
fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
for ax, key, unit in zip(axes.flat, ['v', 'i', 'temperature_C', 'energy'], ['Capacitor voltage (V)', 'Inductor current (A)', 'Resistor temperature (°C)', 'Stored electrical energy (J)']):
    ax.plot(column(positive, 'elapsed_s'), column(positive, 'true_' + key), color='.55', lw=3, label='Synthetic truth')
    for sign, color, label in [('1', '#0072B2', 'Positive fitted branch'), ('-1', '#D55E00', 'Negative fitted branch')]:
        rows = [r for r in paths if r['branch'] == sign]
        ax.plot(column(rows, 'elapsed_s'), column(rows, 'discrete_' + key), color=color, label=label)
        if key != 'energy':
            ax.plot(column(rows, 'elapsed_s'), column(rows, 'analytical_' + key), '--', color=color, lw=1)
    ax.set(xlabel='Elapsed time from t₀ (s)', ylabel=unit)
    if key == 'temperature_C':
        rows = [r for r in read('observations.csv') if r['design'] == '0']
        ax.errorbar(column(rows, 'elapsed_s'), 20 + column(rows, 'value'), yerr=column(rows, 'scale'), fmt='ko', ms=4, label='Temperature observations ± bound')
    ax.legend(fontsize=8)
fig.suptitle('RLC discharge observed through resistor heating | θ=.5, N=128\nDashed: analytical fitted states; two polarities give the same heating and energy')
save(fig, 'reconstruction')

fig, axes = plt.subplots(1, 2, figsize=(13, 6), layout='constrained')
for sign, color, label in [('1', '#0072B2', 'Positive branch'), ('-1', '#D55E00', 'Negative branch')]:
    rows = [r for r in paths if r['branch'] == sign]
    axes[0].plot(column(rows, 'discrete_v'), column(rows, 'discrete_i'), color=color, label=label)
    axes[0].plot(column(rows, 'analytical_v'), column(rows, 'analytical_i'), '--', color=color)
    axes[0].scatter(float(rows[0]['discrete_v']), float(rows[0]['discrete_i']), color=color, marker='o')
    axes[0].scatter(float(rows[-1]['discrete_v']), float(rows[-1]['discrete_i']), color=color, marker='s')
    axes[1].plot(column(rows, 'discrete_energy'), column(rows, 'discrete_temperature_C'), color=color, label=label)
axes[0].set(xlabel='Capacitor voltage (V)', ylabel='Inductor current (A)', title='Electrical phase plane | circle: initial; square: final')
axes[1].set(xlabel='Stored electrical energy (J)', ylabel='Resistor temperature (°C)', title='Thermal versus energy paths overlap')
for ax in axes:
    ax.legend()
save(fig, 'phase')

boxes = read('feasible_boxes.csv')
fig, axes = plt.subplots(2, 3, figsize=(17, 10), layout='constrained')
for j, (case, title) in enumerate(zip(['temperature_only', 'signed_voltage', 'weak_signed_voltage'], ['Temperature only: both polarities', 'Added ±0.02 V datum: one polarity', 'Added ±3 V datum: ambiguity remains'])):
    for row in (0, 1):
        ax = axes[row, j]
        positive_rows = []
        for which, color in [('numerical_outer', '#0072B2'), ('analytical', '#D55E00')]:
            rows = [r for r in boxes if r['case'] == case and r['set'] == which and (row == 0 or float(r['lo_v']) > 0)]
            vertices = [[(float(r['lo_v']), float(r['lo_i'])), (float(r['hi_v']), float(r['lo_i'])), (float(r['hi_v']), float(r['hi_i'])), (float(r['lo_v']), float(r['hi_i']))] for r in rows]
            ax.add_collection(PolyCollection(vertices, facecolors=color, edgecolors='none', alpha=.35, label=which))
            if row == 1:
                positive_rows.extend(rows)
        ax.scatter([1.2], [.7], c='k', marker='*', s=70, label='Truth')
        ax.set(xlabel='Initial capacitor voltage (V)', ylabel='Initial inductor current (A)', title=title if row == 0 else 'Positive region: detail only')
        if row == 0:
            ax.set(xlim=(-1.6, 1.6), ylim=(-1.4, 1.4))
        else:
            ax.set(xlim=(min(float(r['lo_v']) for r in positive_rows) - .005, max(float(r['hi_v']) for r in positive_rows) + .005), ylim=(min(float(r['lo_i']) for r in positive_rows) - .005, max(float(r['hi_i']) for r in positive_rows) + .005))
        ax.legend(fontsize=8)
fig.suptitle('Global bounded feasible covers | entire stated prior searched; ordinary double precision\nTop: full prior; bottom: magnified positive regions (negative regions remain in qualification)')
save(fig, 'feasible_sets')

errors, rates = read('convergence.csv'), read('orders.csv')
keys = [k for k in errors[0] if k not in ('design', 'theta', 'steps', 'h_s')]
fig, axes = plt.subplots(3, 5, figsize=(19, 11), layout='constrained')
for ax, key in zip(axes.flat, keys):
    for design, theta, color, marker, style in [(0, .5, '#0072B2', 'o', '-'), (0, .65, '#D55E00', 's', '-'), (1, .5, '#009E73', '^', '--'), (1, .65, '#CC79A7', 'v', '--')]:
        rows = [r for r in errors if int(r['design']) == design and float(r['theta']) == theta]
        order = next(float(r['order']) for r in rates if int(r['design']) == design and float(r['theta']) == theta and r['quantity'] == key)
        ax.loglog(column(rows, 'h_s'), column(rows, key), marker=marker, linestyle=style, color=color, ms=3, label=f'{"Signed" if design else "Thermal"}, θ={theta}: {order:.3f}')
    ax.set(title=key.replace('_', ' '), xlabel='Step size (s)', ylabel='Error (stated scales)')
    ax.legend(fontsize=6)
for ax in list(axes.flat)[len(keys):]:
    ax.set_visible(False)
fig.suptitle('Independent analytical references | thirteen measured temporal-convergence quantities\nThe temperature-only initial-control comparison follows the positive branch')
save(fig, 'convergence')

bounds, decisions = read('qualification_envelopes.csv'), read('qualification_decisions.csv')
cases = [('temperature_only', 'Temperature-only data'), ('signed_voltage', 'Added signed voltage | ±0.02 V'), ('weak_signed_voltage', 'Added signed voltage | ±3 V')]
goals = ['voltage', 'current', 'temperature_rise', 'stored_energy']
titles = ['Voltage', 'Current', 'Temperature rise', 'Stored energy']
units = ['V', 'A', 'K rise', 'J']
fig, axes = plt.subplots(3, 4, figsize=(17, 11), layout='constrained')
for i, (case, label) in enumerate(cases):
    for j, goal in enumerate(goals):
        ax = axes[i, j]
        rows = [r for r in bounds if r['case'] == case and r['goal'] == goal]
        d = next(r for r in decisions if r['case'] == case and r['goal'] == goal)
        ax.semilogy(column(rows, 'elapsed_s'), column(rows, 'global_bound') / column(rows, 'tolerance'), color='#0072B2', label='All feasible branches + numerical error')
        ax.semilogy(column(rows, 'elapsed_s'), column(rows, 'selected_branch_bound') / column(rows, 'tolerance'), ':', color='#D55E00', label='Selected positive branch only (diagnostic)')
        ax.axhline(1, c='k', lw=1)
        verdict, color = 'NOT QUALIFIED', '#B74425'
        if d['qualified'] == '1':
            time = float(d['persistent_s'])
            verdict, color = f'QUALIFIED from {time:.8g} s', '#14814B'
            ax.axvline(time, color=color, linestyle='--')
        ax.set_title(f'{titles[j]} | tolerance {float(d["tolerance"]):g} {units[j]}\n{verdict}', fontsize=10, color=color)
        if j == 0:
            ax.set_ylabel(label + '\nError bound / tolerance', fontsize=9)
        if i == 2:
            ax.set_xlabel('Elapsed time (s)')
        if i == 0 and j == 0:
            ax.legend(fontsize=7)
fig.suptitle('Engineering qualification | θ=.5, N=1024; whole intervals through 8 elapsed seconds\nDotted curves show the false confidence obtained by discarding the other feasible polarity')
save(fig, 'qualification')

fig, ax = plt.subplots(figsize=(15, 4.5), layout='constrained')
ax.axis('off')
rows, colors = [], []
for case, label in cases:
    ds = [next(r for r in decisions if r['case'] == case and r['goal'] == goal) for goal in goals]
    row = [label] + [f'From {float(r["persistent_s"]):.8g} s' if r['qualified'] == '1' else 'NOT QUALIFIED' for r in ds]
    row.append(f'From {float(ds[0]["joint_s"]):.8g} s' if ds[0]['joint_qualified'] == '1' else 'NOT QUALIFIED')
    rows.append(row)
    colors.append(['#F0F2F3'] + ['#D7EBDF' if x.startswith('From') else '#F4DEDA' for x in row[1:]])
table = ax.table(cellText=rows, cellColours=colors, colLabels=['Observations'] + titles + ['ALL GOALS'], cellLoc='center', loc='center', colWidths=[.28, .144, .144, .144, .144, .144])
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 3)
ax.set_title('ELECTROTHERMAL ENGINEERING QUALIFICATION\nGood local fits do not establish global polarity', fontsize=15, pad=20)
fig.text(.5, .015, 'Retrospective decisions under bounded measurement errors and the full stated initial-state domain.', ha='center')
save(fig, 'qualification_verdict')

rows = read('costates.csv')
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
for ax, key, label in zip(axes, ['v', 'i'], ['Capacitor-voltage', 'Inductor-current']):
    ax.plot(column(rows, 'elapsed_s'), column(rows, 'analytical_' + key), 'k--', lw=2, label='Analytical')
    ax.plot(column(rows, 'elapsed_s'), column(rows, 'discrete_' + key), color='#0072B2', label='Discrete cost-to-go')
    ax.set(xlabel='Elapsed time (s)', ylabel=f'Dissipated-energy sensitivity ({"J/V" if key == "v" else "J/A"})', title=label + ' running-cost costate')
    ax.legend()
save(fig, 'costates')
print('Rendered electrothermal C++ outputs; Python performs no forward solve, optimization, or qualification.')
