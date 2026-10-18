"""Plot C++20 qualification outputs. No inference or qualification decisions."""
from pathlib import Path
import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
out=root/'results'/'qualification'
def read(name,numeric=True):
    with (out/name).open() as f:
        rows=list(csv.DictReader(f))
    if numeric:return [{k:float(v) for k,v in row.items()} for row in rows]
    return rows
curves=read('curves.csv');added=read('added_exposure.csv');refinement=read('refinement.csv')
steps=read('continuous_envelopes.csv',False)
times={(r['case'],r['goal']):r for r in read('qualification_times.csv',False)}
decisions={r['case']:r for r in read('qualification_decisions.csv',False)}
def verdict(case):
    row=decisions[case]
    if row['all_goals_status']=='qualified':
        return f"ALL REQUESTED GOALS QUALIFIED from {float(row['all_goals_elapsed_min']):.5g} min"
    return 'NOT ALL GOALS QUALIFIED (blocking goal: '+row['blocking_goals']+')'
def decision_cell(row,goal):
    if row[goal+'_status']=='qualified':return f"QUALIFIED\nfrom {float(row[goal+'_elapsed_min']):.5g} min"
    return 'NOT QUALIFIED'
def time_label(case,goal):
    row=times[case,goal]
    return f"tQ={float(row['elapsed_min']):.5g} min" if row['status']=='qualified' else 'unqualified on [0, 10] min'
def mark_time(ax,case,goal,color):
    row=times[case,goal]
    if row['status']=='qualified':
        ax.axvline(float(row['elapsed_min']),color=color,ls='-.',lw=1.2,alpha=.8)
plt.rcParams.update({'font.size':10})

fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for ax,component,tolerance,title in zip(axes,['x','q'],[.02,.025],['Concentration','Accumulated exposure']):
    ax.plot([r['elapsed'] for r in curves],[r[f'exact_bound_{component}']/tolerance for r in curves],'k-',label='Analytical: '+time_label('baseline_analytical',component))
    ax.plot([r['elapsed'] for r in curves],[r[f'raw_bound_{component}']/tolerance for r in curves],color='tab:orange',ls=':',label='Raw nodes (uncorrected): '+time_label('baseline_raw_discrete',component))
    for case,color,label in [('baseline','tab:blue','Physical continuous envelope, baseline'),('added','tab:green','Physical continuous envelope, added exposure')]:
        rr=[r for r in steps if r['case']==case]
        starts=[float(r['start']) for r in rr]+[float(rr[-1]['end'])]
        values=[float(r[f'physical_bound_{component}'])/tolerance for r in rr]
        key='baseline_physical_continuous' if case=='baseline' else 'added_physical_continuous'
        ax.stairs(values,starts,color=color,label=label+'\n'+time_label(key,component),baseline=None)
        mark_time(ax,key,component,color)
    mark_time(ax,'baseline_analytical',component,'black')
    mark_time(ax,'baseline_raw_discrete',component,'tab:orange')
    ax.axhline(1,color='red',ls='--',label='Tolerance')
    ax.set(title=title,xlabel='Elapsed time t − t₀ (min)',ylabel='Bounded-error envelope / tolerance')
    ax.legend(fontsize=7);ax.grid(alpha=.2)
fig.suptitle('Physical qualification verdict (C++20, bounded-error data)\nBaseline: '+verdict('baseline')+'\nAdded exposure: '+verdict('added_exposure'),fontsize=10)
fig.savefig(out/'qualification.png',dpi=160);plt.close(fig)

fig,ax=plt.subplots(figsize=(12,4.5),layout='constrained');ax.axis('off')
table_rows=[]
for case,label in [('baseline','Baseline\nconcentration data'),('added_exposure','Added exposure\nrefined grid')]:
    row=decisions[case]
    table_rows.append([label,decision_cell(row,'concentration'),decision_cell(row,'exposure'),decision_cell(row,'all_goals'),row['blocking_goals'] or 'none'])
table=ax.table(cellText=table_rows,colLabels=['Case','Concentration','Exposure','All requested goals','Blocking goal'],cellLoc='center',loc='center',colWidths=[.20,.20,.20,.24,.16])
table.auto_set_font_size(False);table.set_fontsize(11);table.scale(1,3.8)
for (i,j),cell in table.get_celld().items():
    if i==0:cell.set_facecolor('#e8edf2');cell.set_text_props(weight='bold')
    elif j in [1,2,3]:cell.set_facecolor('#f6d6d5' if 'NOT QUALIFIED' in cell.get_text().get_text() else '#dcefdc')
ax.set_title('Physical qualification determination\nNumerical-error-aware envelopes covering the full interval between nodes',fontsize=14,pad=20)
ax.text(.5,.08,'Times are elapsed from t₀; qualification persists through elapsed 10 min.\nAdded-measurement results use retrospective reconstruction with all observations.\nPassing numerical tests validates these decisions; it does not mean every case qualifies.',ha='center',va='center',transform=ax.transAxes,fontsize=10)
fig.savefig(out/'qualification_verdict.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for j,component in enumerate(['x','q']):
    ax=axes[0,j]
    ax.plot([r['elapsed'] for r in curves],[r[f'{component}_nominal'] for r in curves],label='C++ fit, Gaussian record')
    ax.plot([r['elapsed'] for r in curves],[r[f'truth_{component}'] for r in curves],'k--',label='Synthetic truth')
    ax.fill_between([r['elapsed'] for r in curves],[r[f'{component}_nominal']-r[f'local95_{component}'] for r in curves],[r[f'{component}_nominal']+r[f'local95_{component}'] for r in curves],alpha=.2,label='Local linear 95% indicator')
    ax.set(xlabel='Elapsed time (min)',ylabel='Normalized concentration x' if j==0 else 'Exposure q (min)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1,j];tol=.02 if j==0 else .025
    ax.plot([r['elapsed'] for r in curves],[r[f'local95_{component}']/tol for r in curves],label='Local linear 95%: '+time_label('local_gaussian',component))
    ax.plot([r['elapsed'] for r in curves],[r[f'joint95_physical_bound_{component}']/tol for r in curves],label='Joint-noise-box 95%: '+time_label('joint95_gaussian',component))
    mark_time(ax,'local_gaussian',component,'tab:blue')
    mark_time(ax,'joint95_gaussian',component,'tab:orange')
    ax.axhline(1,color='red',ls='--');ax.set(xlabel='Elapsed time (min)',ylabel='Statistical allowance / tolerance');ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.savefig(out/'information.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for ax,key,title in zip(axes,['control_error','goal_error','hessian_error'],['Reconstructed initial state','Exposure goal adjoint','Reduced Hessian']):
    for theta in [.5,.65]:
        rr=[r for r in refinement if abs(r['theta']-theta)<1e-12]
        ax.loglog([r['h'] for r in rr],[r[key] for r in rr],'o-',label=f'θ={theta:g}')
    ax.set(title=title,xlabel='Maximum step h (min)',ylabel='Absolute error to analytical reference');ax.legend();ax.grid(which='both',alpha=.2)
fig.savefig(out/'convergence.png',dpi=160);plt.close(fig)

fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
for theta,color in [(.5,'tab:blue'),(.65,'tab:orange')]:
    rr=[r for r in refinement if abs(r['theta']-theta)<1e-12 and math.isfinite(r['physical_qualification_x'])]
    ax.semilogx([r['h'] for r in rr],[r['physical_qualification_x'] for r in rr],'o-',color=color,label=f'Physical node qualification, θ={theta:g}')
reference=float(times['baseline_analytical','x']['elapsed_min'])
ax.axhline(reference,color='black',ls='--',label=f'Analytical continuous tQ={reference:.5f} min')
ax.set(xlabel='Maximum step h (min)',ylabel='Persistent concentration qualification time (elapsed min)',title='Qualification time under mesh refinement')
ax.legend(fontsize=8);ax.grid(which='both',alpha=.2)
ax.text(.98,.03,'Unqualified coarse meshes have no finite tQ and are omitted',transform=ax.transAxes,ha='right',fontsize=8)
fig.savefig(out/'qualification_time_convergence.png',dpi=160);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for j,component in enumerate(['x','q']):
    axes[j].plot([r['elapsed'] for r in curves],[r[f'raw_bound_{component}'] for r in curves],color='k',label='Exact nonlinear discrete interval')
    axes[j].plot([r['elapsed'] for r in curves],[r[f'linear_bound_{component}'] for r in curves],ls=':',label='Linearized box support')
    axes[j].plot([r['elapsed'] for r in curves],[r[f'taylor_bound_{component}'] for r in curves],ls='--',label='Linear plus curvature remainder')
    axes[j].set(xlabel='Elapsed time (min)',ylabel='Concentration uncertainty' if j==0 else 'Exposure uncertainty (min)',title='Nonlinear concentration envelope' if j==0 else 'Nonlinear exposure envelope')
    axes[j].legend(fontsize=8);axes[j].grid(alpha=.2)
fig.savefig(out/'nonlinear_remainder.png',dpi=160);plt.close(fig)
print('Plotted C++ qualification outputs; no numerical model or qualification logic runs in Python.')
