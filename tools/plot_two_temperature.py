# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Render thermal results produced entirely by C++20."""
from pathlib import Path
import csv,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1];out=root/'results/two_temperature'
def read(name,numeric=True):
    with (out/name).open() as f:rows=list(csv.DictReader(f))
    return [{k:float(v) for k,v in row.items()} for row in rows] if numeric else rows
curves=read('curves.csv');observations=read('observations.csv');convergence=read('convergence.csv')
steps=read('continuous_envelopes.csv',False);decisions={r['case']:r for r in read('qualification_decisions.csv',False)}
summary={r['quantity']:float(r['value']) for r in read('summary.csv',False)}
plt.rcParams.update({'font.size':10})

fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for ax,key,exact,title in zip(axes,['fluid_C','wall_C'],['truth_fluid_C','truth_wall_C'],['Fluid temperature','Wall temperature']):
    ax.plot([r['elapsed'] for r in curves],[r[key] for r in curves],label='C++ adjoint-reconstructed trajectory')
    ax.plot([r['elapsed'] for r in curves],[r[exact] for r in curves],'k--',label='Analytical synthetic truth')
    if key=='fluid_C':ax.errorbar([r['elapsed'] for r in observations],[r['observation_fluid_C'] for r in observations],yerr=.1,fmt='o',ms=3,label='Fluid observations ±1 SD')
    ax.set(title=title,xlabel='Elapsed time t − t₀ (min)',ylabel='Temperature (°C)');ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.savefig(out/'reconstruction.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
labels={'primary':'Primary coarse grid','refined':'Refined grid','added_wall_sensor':'Refined + wall sensor'}
colors={'primary':'tab:orange','refined':'tab:blue','added_wall_sensor':'tab:green'}
for j,(goal,tolerance,title) in enumerate([('fluid',.2,'Fluid-temperature error'),('energy',60.,'Stored-energy error')]):
    ax=axes[j]
    for case in labels:
        rows=[r for r in steps if r['case']==case];edges=[float(r['start']) for r in rows]+[float(rows[-1]['end'])]
        value_key='fluid_bound_K' if goal=='fluid' else 'energy_bound_kJ'
        values=[float(r[value_key])/tolerance for r in rows]
        row=decisions[case];qualified=row[goal+'_status']=='qualified'
        time=float(row[goal+'_elapsed_min']) if qualified else None
        caption=f'tQ={time:.5g} min' if qualified else 'NOT QUALIFIED'
        ax.stairs(values,edges,color=colors[case],baseline=None,label=labels[case]+': '+caption)
        if qualified:ax.axvline(time,color=colors[case],ls='-.',lw=1)
    ax.plot([r['elapsed'] for r in curves],[r['exact_ellipse_fluid_K' if j==0 else 'exact_ellipse_energy_kJ']/tolerance for r in curves],'k:',label='Analytical 95% ellipse reference (different set)')
    ax.axhline(1,color='red',ls='--',label='Tolerance');ax.set_yscale('log')
    ax.set(title=title,xlabel='Elapsed time t − t₀ (min)',ylabel='Whole-step error allowance / tolerance');ax.grid(which='both',alpha=.2);ax.legend(fontsize=7)
fig.suptitle('95% noise-box physical qualification: numerical error included; initial temperatures in [0,60] °C',fontsize=11)
fig.savefig(out/'qualification.png',dpi=160);plt.close(fig)

fig,ax=plt.subplots(figsize=(12,5),layout='constrained');ax.axis('off')
def cell(row,goal):return 'QUALIFIED\nfrom '+f"{float(row[goal+'_elapsed_min']):.6g}"+' min' if row[goal+'_status']=='qualified' else 'NOT QUALIFIED'
data=[]
for case in labels:
    row=decisions[case];data.append([labels[case],cell(row,'fluid'),cell(row,'energy'),cell(row,'all_goals'),row['blocking_goals'] or 'none'])
table=ax.table(cellText=data,colLabels=['Case','Fluid temperature','Stored energy','All requested goals','Blocking goals'],colWidths=[.2,.2,.2,.24,.16],cellLoc='center',loc='center')
table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,3.2)
for (i,j),c in table.get_celld().items():
    if i==0:c.set_facecolor('#e8edf2');c.set_text_props(weight='bold')
    elif j in [1,2,3]:c.set_facecolor('#f6d6d5' if 'NOT QUALIFIED' in c.get_text().get_text() else '#dcefdc')
ax.set_title('Two-temperature physical qualification verdict\nFluid tolerance 0.2 K; stored-energy tolerance 60 kJ',fontsize=14,pad=15)
ax.text(.5,.03,'Times are elapsed from t₀ and qualification persists through elapsed 10 min.\n95% simultaneous independent-Gaussian noise box, known initial domain [0,60] °C for both temperatures.\nWhole-step numerical error included; ordinary floating-point bounds. Tests passing does not imply every case qualifies.',ha='center',transform=ax.transAxes,fontsize=9)
fig.savefig(out/'qualification_verdict.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(2,4,figsize=(14,7),layout='constrained')
keys=['state_error','control_error','tangent_error','goal_error','covariance_error','running_cost_error','running_gradient_error']
for ax,key in zip(axes.flat,keys):
    for theta in [.5,.65]:
        rr=[r for r in convergence if abs(r['theta']-theta)<1e-12]
        order=math.log2(rr[-2][key]/rr[-1][key]);ax.loglog([r['h'] for r in rr],[r[key] for r in rr],'o-',label=f'θ={theta:g}, order={order:.4f}')
    ax.set(title=key.replace('_',' ').capitalize(),xlabel='Maximum step h (min)',ylabel='Absolute analytical-reference error');ax.legend(fontsize=8);ax.grid(which='both',alpha=.2)
axes.flat[-1].axis('off');axes.flat[-1].text(.05,.8,'C++20 theta optimality library\n\nTwo initial-temperature controls\nExplicit nodal observation costs\nFull discrete tangent + adjoint\nRunning theta cost verification\nGoal information + covariance\nConservative energy balance',va='top')
fig.savefig(out/'convergence.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for j,goal in enumerate(['fluid','energy']):
    exact_key='exact_ellipse_fluid_K' if j==0 else 'exact_ellipse_energy_kJ';discrete_key='discrete_ellipse_fluid_K' if j==0 else 'discrete_ellipse_energy_kJ';tol=.2 if j==0 else 60.
    axes[j].plot([r['elapsed'] for r in curves],[r[exact_key]/tol for r in curves],label='Analytical simultaneous 95% ellipse')
    axes[j].plot([r['elapsed'] for r in curves],[r[discrete_key]/tol for r in curves],'--',label='Discrete-model ellipse (bias excluded)')
    tq=summary['analytical_ellipse_'+goal+'_qualification_min'];axes[j].axvline(tq,color='k',ls=':',label=f'Analytical tQ={tq:.5g} min')
    axes[j].axhline(1,color='red',ls='--');axes[j].set(title='Fluid temperature' if j==0 else 'Stored energy',xlabel='Elapsed time (min)',ylabel='Ellipse allowance / tolerance');axes[j].grid(alpha=.2);axes[j].legend(fontsize=7)
fig.suptitle('Retained Gaussian-ellipse reference: distinct from the physical noise-box verdict',fontsize=10)
fig.savefig(out/'ellipse_reference.png',dpi=160);plt.close(fig)
print('Rendered two-temperature C++ outputs; no fitting or qualification decisions run in Python.')
