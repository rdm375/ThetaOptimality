#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Render C++ oscillator outputs. No state, fitting, or convergence logic."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]/'results'/'oscillator'
def rows(name):
    with (root/name).open() as f:return list(csv.DictReader(f))
def col(records,key):return np.array([float(r[key]) for r in records])
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
r=rows('trajectories.csv');obs=rows('observations.csv');t=col(r,'elapsed_s')
fig,axes=plt.subplots(2,2,figsize=(12,8))
for ax,key,label in [(axes[0,0],'q_m','Displacement (m)'),(axes[0,1],'v_m_s','Velocity (m/s)')]:
    ax.plot(t,col(r,'true_'+key),color='0.6',lw=2,label='Synthetic truth')
    ax.plot(t,col(r,'analytical_fit_'+key),'k--',label='Analytical fit')
    ax.plot(t,col(r,'discrete_'+key),color='#2876a5',label='Discrete fit, θ=0.5, N=128')
    ax.set_xlabel('Elapsed time from t₀ (s)');ax.set_ylabel(label);ax.grid(alpha=.2)
axes[0,0].errorbar(col(obs,'elapsed_s'),col(obs,'observed_displacement_m'),yerr=col(obs,'sigma_m'),fmt='o',ms=4,color='#a6462d',capsize=3,label='Displacement data ± scale')
axes[0,0].legend(fontsize=8);axes[0,1].set_title('Hidden velocity: reconstructed without a velocity sensor')
axes[1,0].plot(t,1e3*(col(r,'discrete_q_m')-col(r,'analytical_fit_q_m')),label='Displacement difference (mm)')
velocity_axis=axes[1,0].twinx()
velocity_axis.plot(t,1e3*(col(r,'discrete_v_m_s')-col(r,'analytical_fit_v_m_s')),color='#ba5c2a',label='Velocity difference (mm/s)')
velocity_axis.set_ylabel('Velocity: discrete fit − analytical fit (mm/s)',color='#ba5c2a')
axes[1,0].set_ylabel('Displacement: discrete fit − analytical fit (mm)');axes[1,0].set_xlabel('Elapsed time from t₀ (s)');axes[1,0].legend(fontsize=8);axes[1,0].grid(alpha=.2)
axes[1,1].plot(t,col(r,'true_energy_J'),color='0.6',label='Truth')
axes[1,1].plot(t,col(r,'discrete_energy_J'),color='#2876a5',label='Discrete fitted state')
axes[1,1].set_xlabel('Elapsed time from t₀ (s)');axes[1,1].set_ylabel('Mechanical energy (J)');axes[1,1].legend();axes[1,1].grid(alpha=.2)
fig.suptitle('Displacement-only state finding | known m, c, k | unknown initial displacement and velocity')
fig.tight_layout(rect=(0,0,1,.95));fig.savefig(root/'reconstruction.png',dpi=170);plt.close(fig)
conv=rows('convergence.csv');orders={(float(r['theta']),r['component']):float(r['order']) for r in rows('orders.csv')}
keys=['state_q','state_v','tangent','nodal_adjoint','control_q','control_v','nodal_cost','nodal_gradient','nodal_hessian','running_cost','running_gradient','running_costate','running_hessian']
titles=['Displacement state','Velocity state','Tangent propagation','Observation adjoint','Fitted initial displacement','Fitted initial velocity','Observation objective','Observation gradient','Observation Hessian','Running-cost quadrature','Running-cost gradient','Running-cost costate','Running-cost Hessian']
fig,axes=plt.subplots(4,4,figsize=(15,12))
for ax,key,title in zip(axes.flat,keys,titles):
    for theta,color,marker in [('0.5','#2876a5','o'),('0.65','#ba5c2a','s')]:
        subset=[r for r in conv if float(r['theta'])==float(theta)];h=col(subset,'h_s');err=col(subset,key);order=orders[(float(theta),key)]
        ax.loglog(h,err,marker+'-',color=color,label=f'θ={theta}; order {order:.3f}',ms=4)
    ax.set_title(title,fontsize=10);ax.set_xlabel('Step size h (s)');ax.set_ylabel('Normalized error');ax.grid(which='both',alpha=.2);ax.legend(fontsize=7)
for ax in list(axes.flat)[len(keys):]:ax.set_visible(False)
fig.suptitle('Temporal convergence against analytical references | orders from final two grids',fontsize=15)
fig.tight_layout(rect=(0,0,1,.96));fig.savefig(root/'convergence.png',dpi=160);plt.close(fig)
controls=rows('control_convergence.csv');fig,axes=plt.subplots(1,2,figsize=(11,4.5))
for ax,key,ref,label in [(axes[0],'discrete_q0_m','analytical_q0_m','Initial displacement (m)'),(axes[1],'discrete_v0_m_s','analytical_v0_m_s','Initial velocity (m/s)')]:
    for theta in ['0.5','0.65']:
        subset=[r for r in controls if float(r['theta'])==float(theta)];ax.plot(col(subset,'h_s'),col(subset,key),'o-',label='Discrete θ='+theta)
    ax.axhline(float(controls[0][ref]),color='k',ls='--',label='Analytical noisy-data fit')
    ax.set_xlabel('Step size h (s)');ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Recovered controls converge to the analytical fit of the SAME observations')
fig.tight_layout(rect=(0,0,1,.92));fig.savefig(root/'control_convergence.png',dpi=170);plt.close(fig)
blind=rows('hidden_mode.csv');design=rows('measurement_design.csv');fig,axes=plt.subplots(1,2,figsize=(12,4.8))
for key,label in [('q_vminus_m','Initial v = −0.03 m/s'),('q_vtrue_m','Initial v = 0.12 m/s'),('q_vplus_m','Initial v = 0.25 m/s')]:axes[0].plot(col(blind,'elapsed_s'),col(blind,key),label=label)
poor=[r for r in design if r['design']=='half_period'];tp=col(poor,'elapsed_s');axes[0].plot(tp,col(poor,'true_displacement_m'),'ko',label='Identical half-period samples')
axes[0].set_xlabel('Elapsed time from t₀ (s)');axes[0].set_ylabel('Displacement (m)');axes[0].legend(fontsize=8);axes[0].set_title('Different hidden velocities; identical sampled displacement')
for label,color in [('half_period','#a6462d'),('well_spaced','#2876a5')]:
    subset=[r for r in design if r['design']==label];axes[1].scatter(.02*col(subset,'q0_sensitivity'),.12*col(subset,'v0_sensitivity'),label=label.replace('_',' '),color=color)
axes[1].axhline(0,color='0.7',lw=.7);axes[1].axvline(0,color='0.7',lw=.7)
axes[1].set_xlabel('Sensitivity to q₀ / 0.02 m (m)');axes[1].set_ylabel('Sensitivity to v₀ / 0.12 m/s (m)');axes[1].set_title('Scaled sensitivity rows\nHalf-period rows lie on one line');axes[1].legend()
for ax in axes:ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(root/'measurement_design.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4.5))
for ax,key,label in [(axes[0],'q','λq (m·s)'),(axes[1],'v','λv (m·s²)')]:
    ax.plot(t,col(r,'exact_probe_running_costate_'+key),'k--',label='Analytical costate')
    ax.plot(t,col(r,'discrete_probe_running_costate_'+key),color='#2876a5',label='Corrected discrete costate')
    ax.set_xlabel('Elapsed time from t₀ (s)');ax.set_ylabel(label);ax.grid(alpha=.2);ax.legend()
fig.suptitle('Separate verification functional Jrun = ½ ∫ q² dt | terminal costate zero')
fig.tight_layout(rect=(0,0,1,.92));fig.savefig(root/'costate.png',dpi=170);plt.close(fig)

# Phase portraits use the same C++ analytical and optimized discrete trajectories.
fig,axes=plt.subplots(1,2,figsize=(13,5.7))
ax=axes[0]
ax.plot(col(r,'true_q_m'),col(r,'true_v_m_s'),color='0.65',lw=2.5,label='Analytical synthetic truth')
ax.plot(col(r,'analytical_fit_q_m'),col(r,'analytical_fit_v_m_s'),'k--',lw=1.8,label='Analytical optimal fit (same noisy data)',zorder=4)
ax.plot(col(r,'discrete_q_m'),col(r,'discrete_v_m_s'),color='#2876a5',lw=1.3,marker='.',markevery=8,label='Discrete optimal reconstruction, N=128')
for qkey,vkey,color in [('analytical_fit_q_m','analytical_fit_v_m_s','k'),('discrete_q_m','discrete_v_m_s','#2876a5')]:
    ax.plot(float(r[0][qkey]),float(r[0][vkey]),'o',mfc='white',mec=color,ms=7,zorder=6)
    ax.plot(float(r[-1][qkey]),float(r[-1][vkey]),'s',mfc=color,mec=color,ms=5,zorder=6)
for j in [10,36,65,96]:
    ax.annotate('',xy=(float(r[j+2]['discrete_q_m']),float(r[j+2]['discrete_v_m_s'])),xytext=(float(r[j]['discrete_q_m']),float(r[j]['discrete_v_m_s'])),arrowprops={'arrowstyle':'->','color':'#2876a5','lw':1.5})
ax.plot(0,0,'+',color='0.3',ms=9)
ax.set_xlabel('Displacement q (m)');ax.set_ylabel('Velocity v (m/s)')
ax.set_title('Physical phase plane: damped motion toward equilibrium\nOpen circles: start; squares: end; arrows: increasing time',fontsize=10)
ax.legend(fontsize=8,loc='lower left');ax.grid(alpha=.2)
ax=axes[1]
dq=1e3*(col(r,'discrete_q_m')-col(r,'analytical_fit_q_m'))
dv=1e3*(col(r,'discrete_v_m_s')-col(r,'analytical_fit_v_m_s'))
ax.plot(dq,dv,color='0.7',lw=1,zorder=1)
points=ax.scatter(dq,dv,c=t,cmap='viridis',s=13,zorder=2)
ax.plot(dq[0],dv[0],'o',mfc='white',mec='k',ms=8,label='Initial-state difference',zorder=4)
ax.plot(dq[-1],dv[-1],'ks',ms=6,label='Final-state difference',zorder=4)
ax.axhline(0,color='0.6',lw=.7);ax.axvline(0,color='0.6',lw=.7)
ax.set_xlabel('Discrete − analytical optimal displacement (mm)')
ax.set_ylabel('Discrete − analytical optimal velocity (mm/s)')
ax.set_title('Difference phase plane: reveals the small mismatch\nIncludes both optimized initial-state and propagation errors',fontsize=10)
ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.colorbar(points,ax=ax,label='Elapsed time from t₀ (s)',pad=.02)
phase_status=next(x for x in rows('qualification_decisions.csv') if x['case']=='bounded_primary')['joint_qualified']
fig.suptitle('Analytical solution and optimal discrete reconstruction in phase space | θ=0.5',fontsize=13)
fig.text(.5,.01,'Engineering status of this N=128 case: '+('QUALIFIED' if phase_status=='1' else 'NOT QUALIFIED')+' under the stated bounded-error budget; see the qualification figures.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.035,1,.93));fig.savefig(root/'phase_portraits.png',dpi=190);plt.close(fig)


# All qualification flags and times are computed by C++ and only displayed here.
qd=rows('qualification_decisions.csv');qe=rows('qualification_envelopes.csv')
case_names=['bounded_primary','bounded_refined','bounded_first_order','bounded_first_order_refined','bounded_first_order_4096','bounded_first_order_8192','joint95_refined']
case_labels={'bounded_primary':'Bounded data | θ=0.5, N=128','bounded_refined':'Bounded data | θ=0.5, N=1024','bounded_first_order':'Bounded data | θ=0.65, N=128','bounded_first_order_refined':'Bounded data | θ=0.65, N=1024','bounded_first_order_4096':'Bounded data | θ=0.65, N=4096','bounded_first_order_8192':'Bounded data | θ=0.65, N=8192','joint95_refined':'Joint 95% Gaussian noise box | θ=0.5, N=1024'}
goals=['displacement','velocity','support_force'];goal_labels=['Displacement | tolerance 1 mm','Velocity | tolerance 2 mm/s','Support force | tolerance 0.04 N']
fig,axes=plt.subplots(len(case_names),3,figsize=(15,3.25*len(case_names)),sharex=True)
for i,case in enumerate(case_names):
    for j,goal in enumerate(goals):
        ax=axes[i,j];records=[x for x in qe if x['case']==case and x['goal']==goal];decision=next(x for x in qd if x['case']==case and x['goal']==goal);tol=float(decision['tolerance'])
        intervals=[x for x in records if x['whole_interval_bound']!='']
        ax.step(np.r_[col(intervals,'elapsed_s'),float(intervals[-1]['interval_end_s'])],np.r_[col(intervals,'whole_interval_bound')/tol,float(intervals[-1]['whole_interval_bound'])/tol],where='post',color='#b45b27',label='Physical whole-interval bound')
        ax.plot(col(records,'elapsed_s'),col(records,'analytical_bound')/tol,'--',color='#2876a5',label='Analytical physical reference')
        ax.plot(col(records,'elapsed_s'),col(records,'local95_indicator')/tol,':',color='0.55',label='Gaussian ellipse indicator only')
        ax.axhline(1,color='k',lw=.9)
        if decision['first_interval_time_s']!='' and decision['first_interval_time_s']!=decision['continuous_time_s']:
            ax.axvline(float(decision['first_interval_time_s']),color='#a44b3b',ls=':',alpha=.7,label='First below: not persistent')
        if decision['qualified']=='1':
            time=float(decision['continuous_time_s']);ax.axvline(time,color='#22713d',ls='--',lw=1);status=f'QUALIFIED from {time:.6g} s'
        else:status='NOT QUALIFIED'
        ax.set_yscale('log');ax.set_ylim(.06,max(2.,1.4*np.max(col(intervals,'whole_interval_bound')/tol)));ax.grid(which='both',alpha=.15)
        ax.set_title(goal_labels[j]+'\n'+status,fontsize=9,color='#22713d' if decision['qualified']=='1' else '#9e3226')
        if j==0:ax.set_ylabel(case_labels[case]+'\nError bound / tolerance',fontsize=9)
        if i==len(case_names)-1:ax.set_xlabel('Elapsed time from t₀ (s)')
axes[0,0].legend(fontsize=7,loc='upper right');axes[-1,0].legend(fontsize=7,loc='upper right')
fig.suptitle('Persistent oscillator state qualification | whole intervals through 8 s | numerical error included',fontsize=14)
fig.text(.5,.012,'All-data retrospective qualification. Deterministic bounds and Gaussian noise boxes are different uncertainty assumptions.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.025,1,.96));fig.savefig(root/'qualification.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(13,6.5));ax.axis('off')
columns=['Experiment','Displacement','Velocity','Support force','ALL GOALS'];table_rows=[];colors=[]
for case in case_names:
    ds=[next(x for x in qd if x['case']==case and x['goal']==goal) for goal in goals]
    text=[case_labels[case]];rowcolors=['#f1f3f5']
    for d in ds:
        text.append('From '+d['continuous_time_s']+' s' if d['qualified']=='1' else 'NOT QUALIFIED');rowcolors.append('#d8ecdf' if d['qualified']=='1' else '#f4deda')
    d=ds[0];text.append('From '+d['joint_time_s']+' s' if d['joint_qualified']=='1' else 'NOT QUALIFIED');rowcolors.append('#d8ecdf' if d['joint_qualified']=='1' else '#f4deda');table_rows.append(text);colors.append(rowcolors)
table=ax.table(cellText=table_rows,colLabels=columns,cellColours=colors,colWidths=[.33,.17,.17,.17,.16],loc='center',cellLoc='center');table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,2.5)
ax.set_title('ENGINEERING QUALIFICATION DETERMINATION\n1 mm displacement | 2 mm/s velocity | 0.04 N support-force error',pad=18,fontsize=13)
fig.text(.5,.05,'Times are elapsed from t₀; every remaining whole interval must satisfy the tolerance. Test success is a separate result.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.08,1,1));fig.savefig(root/'qualification_verdict.png',dpi=180);plt.close(fig)

print('Plotted C++ oscillator outputs; Python performs no simulation, fitting, or order estimation.')

# The poor-timing experiment has its own data files and its own conclusion.
# The C++ executable performs fitting, feasible-set propagation, and qualification.
poor_root=Path(__file__).resolve().parents[1]/'results'/'oscillator_poor_timing'
def poor_rows(name):
    with (poor_root/name).open() as f:return list(csv.DictReader(f))
pr=poor_rows('reconstruction.csv');pt=col(pr,'elapsed_s')
fig,axes=plt.subplots(1,3,figsize=(15,4.8))
for ax,truth_key,fit_key,analytic_key,bound_key,label,unit in [
    (axes[0],'true_q_m','discrete_q_m','qualification_center_q_m','displacement_bound_m','Displacement','m'),
    (axes[1],'true_v_m_s','discrete_v_m_s','qualification_center_v_m_s','velocity_bound_m_s','Velocity','m/s')]:
    truth_values=col(pr,truth_key);fit_values=col(pr,fit_key);center=col(pr,analytic_key);bound=col(pr,bound_key)
    ax.plot(pt,truth_values,color='black',lw=1.8,label='Synthetic truth')
    ax.plot(pt,center,'--',color='#2876a5',lw=1.5,label='Theta trajectory from feasible representative')
    ax.plot(pt,fit_values,color='#b45b27',lw=1.1,label='Theta discrete reconstruction')
    ax.fill_between(pt,fit_values-bound,fit_values+bound,color='#b45b27',alpha=.17,label='Bounded feasible-state envelope')
    ax.set_xlabel('Elapsed time from $t_0$ (s)');ax.set_ylabel(f'{label} ({unit})');ax.grid(alpha=.2)
    ax.set_title(label+' reconstruction');ax.legend(fontsize=7)
phase=poor_rows('feasible_sets.csv')
for name,color,label in [('analytical','#2876a5','Exact physical feasible set'),('numerical_outer','#b45b27','Numerically enlarged outer set')]:
    verts=[r for r in phase if r['set']==name]
    q=np.array([float(x['q0_m']) for x in verts]);v=np.array([float(x['v0_m_s']) for x in verts])
    axes[2].fill(np.r_[q,q[0]],np.r_[v,v[0]],alpha=.2,color=color,label=label)
    axes[2].plot(np.r_[q,q[0]],np.r_[v,v[0]],color=color)
axes[2].plot([.02],[.12],'ko',ms=4,label='Synthetic truth')
axes[2].set_xlabel('Initial displacement $q_0$ (m)');axes[2].set_ylabel('Initial velocity $v_0$ (m/s)')
axes[2].set_title('Retained initial-state uncertainty');axes[2].grid(alpha=.2);axes[2].legend(fontsize=7)
qdec=poor_rows('qualification.csv')
labels={'displacement':'q','velocity':'v','support_force':'force'}
for j,d in enumerate(qdec):
    if j==2:
        axes[2].text(.03,.95,'Support force: NOT QUALIFIED',transform=axes[2].transAxes,va='top',fontsize=8,color='#9e3226',bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
        continue
    msg=f"{labels[d['goal']]}: NOT QUALIFIED" if d['qualified']=='0' else f"{labels[d['goal']]}: qualified from {float(d['persistent_time_s']):.3g} s"
    axes[j].text(.03,.04,msg,transform=axes[j].transAxes,fontsize=8,color='#9e3226' if d['qualified']=='0' else '#22713d',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
fig.suptitle('Half-period-only displacement sampling | theta = 0.5, N = 128\nExact sensor map is blind to initial velocity; qualification retains that direction',fontsize=13)
fig.tight_layout(rect=(0,0,1,.91));fig.savefig(poor_root/'poor_timing_reconstruction.png',dpi=180);plt.close(fig)

convp=poor_rows('convergence.csv')
fig,axes=plt.subplots(1,2,figsize=(10,4))
axes[0].loglog(col(convp,'h_s'),col(convp,'fixed_state_error'),'o-',color='#2876a5')
axes[0].invert_xaxis();axes[0].set_xlabel('Maximum nominal step $h$ (s)');axes[0].set_ylabel('Scaled fixed-state error')
axes[0].set_title('Midpoint forward-state convergence');axes[0].grid(which='both',alpha=.2)
axes[1].loglog(col(convp,'h_s'),np.abs(col(convp,'v0_fit_m_s')),'o-',color='#b45b27',label='Magnitude of fitted $v_0$')
axes[1].loglog(col(convp,'h_s'),col(convp,'scaled_information_small_eigenvalue'),'s--',color='#2876a5',label='Small information eigenvalue')
axes[1].invert_xaxis();axes[1].set_xlabel('Maximum nominal step $h$ (s)');axes[1].set_ylabel('Magnitude / scaled information')
axes[1].set_title('Discrete phase leakage under refinement');axes[1].grid(which='both',alpha=.2);axes[1].legend(fontsize=8)
fig.suptitle('The state solver converges; the velocity estimate does not become identifiable')
fig.tight_layout();fig.savefig(poor_root/'poor_timing_convergence.png',dpi=180);plt.close(fig)
print('Plotted the poorly timed oscillator reconstruction and refinement diagnostics.')

# Direct matched comparison: same oscillator, four data, noise realization,
# bounded error, prior, theta method, and qualification tolerances. Only times change.
proper_root=Path(__file__).resolve().parents[1]/'results'/'oscillator_proper_timing'
def proper_rows(name):
    with (proper_root/name).open() as f:return list(csv.DictReader(f))

fig,axes=plt.subplots(3,3,figsize=(15,12))
designs=[
    (poor_root,poor_rows,'Half-period only','#b45b27'),
    (proper_root,proper_rows,'Properly spaced','#2876a5'),
]
for root_dir,reader,label,color in designs:
    reconstruction=reader('reconstruction.csv')
    tau=col(reconstruction,'elapsed_s')
    observations=reader('observations.csv')
    for row,truth_key,disc_key,center_key,bound_key,ylabel in [
        (0,'true_q_m','discrete_q_m','qualification_center_q_m','displacement_bound_m','Displacement (m)'),
        (1,'true_v_m_s','discrete_v_m_s','qualification_center_v_m_s','velocity_bound_m_s','Velocity (m/s)')]:
        ax=axes[row,0 if root_dir==poor_root else 1]
        truth_values=col(reconstruction,truth_key)
        discrete_values=col(reconstruction,disc_key)
        center=col(reconstruction,center_key)
        bound=col(reconstruction,bound_key)
        ax.plot(tau,truth_values,'k-',lw=1.7,label='Analytical truth')
        ax.plot(tau,center,'--',color='#2876a5',lw=1.3,label='Feasible representative')
        ax.plot(tau,discrete_values,color=color,lw=1.1,label='Discrete least-squares fit')
        ax.fill_between(tau,center-bound,center+bound,color=color,alpha=.18,label='Physical feasible envelope')
        if row==0:
            ax.errorbar(col(observations,'elapsed_s'),col(observations,'observed_m'),
                        yerr=col(observations,'measurement_bound_m'),fmt='o',ms=3,
                        color=color,capsize=2,label='Samples ± bound')
        ax.set_title(f'{label}: {"displacement" if row==0 else "velocity"}')
        ax.set_xlabel('Elapsed time from $t_0$ (s)');ax.set_ylabel(ylabel);ax.grid(alpha=.2)
        if row==0 and root_dir==poor_root:ax.legend(fontsize=7,loc='best')

    phase_records=reader('feasible_sets.csv')
    phase_ax=axes[2,0 if root_dir==poor_root else 1]
    for set_name,linestyle,alpha,set_label in [
        ('analytical','-',.18,'Exact physical feasible set'),
        ('numerical_outer','--',.10,'Numerically enlarged outer set')]:
        verts=[r for r in phase_records if r['set']==set_name]
        q=np.array([float(r['q0_m']) for r in verts]);v=np.array([float(r['v0_m_s']) for r in verts])
        phase_ax.fill(np.r_[q,q[0]],np.r_[v,v[0]],color=color,alpha=alpha)
        phase_ax.plot(np.r_[q,q[0]],np.r_[v,v[0]],color=color,ls=linestyle,lw=1,
                      label=f'{label}: {set_label}')
    phase_ax.plot([.02],[.12],'ko',ms=4,label='Synthetic truth')
    phase_ax.set_xlabel('Initial displacement $q_0$ (m)');phase_ax.set_ylabel('Initial velocity $v_0$ (m/s)')
    phase_ax.set_title(f'{label}: retained initial states');phase_ax.grid(alpha=.2)
    phase_ax.legend(fontsize=7)

# Show rank loss under mesh refinement and distinguish it from forward-state convergence.
for root_dir,reader,label,color in designs:
    records=reader('convergence.csv')
    axes[0,2].loglog(col(records,'h_s'),col(records,'scaled_information_small_eigenvalue'),
                     'o-',color=color,label=label)
axes[0,2].invert_xaxis();axes[0,2].set_xlabel('Nominal step $h$ (s)')
axes[0,2].set_ylabel('Small scaled information eigenvalue')
axes[0,2].set_title('Information under mesh refinement');axes[0,2].grid(which='both',alpha=.2)
axes[0,2].legend(fontsize=8)

proper_convergence=proper_rows('convergence.csv');poor_convergence=poor_rows('convergence.csv')
axes[1,2].loglog(col(proper_convergence,'h_s'),col(proper_convergence,'scaled_control_error'),
                 'o-',color='#2876a5',label='Proper: recovered-control error')
axes[1,2].loglog(col(proper_convergence,'h_s'),col(proper_convergence,'fixed_state_error'),
                 's--',color='0.4',label='Proper: fixed-state error')
axes[1,2].loglog(col(poor_convergence,'h_s'),np.abs(col(poor_convergence,'v0_fit_m_s')),
                 '^-',color='#b45b27',label='Poor: $|$fitted $v_0|$')
axes[1,2].invert_xaxis();axes[1,2].set_xlabel('Nominal step $h$ (s)')
axes[1,2].set_ylabel('Scaled error / fitted velocity magnitude')
axes[1,2].set_title('Recovery vs. phase-error leakage');axes[1,2].grid(which='both',alpha=.2)
axes[1,2].legend(fontsize=7)

axes[2,2].axis('off');decision_text=[]
for root_dir,reader,label,_ in designs:
    qd=reader('qualification.csv')
    joint=all(r['qualified']=='1' for r in qd)
    times_qual=[float(r['persistent_time_s']) for r in qd if r['persistent_time_s']]
    status=f'JOINTLY QUALIFIED from {max(times_qual):.3g} s' if joint and times_qual else 'NOT QUALIFIED'
    decision_text.append(f'{label}: {status}')
axes[2,2].text(.04,.88,'Engineering qualification',transform=axes[2,2].transAxes,
               fontsize=13,weight='bold',va='top')
axes[2,2].text(.04,.68,'\n\n'.join(decision_text),transform=axes[2,2].transAxes,
               fontsize=11,va='top',color='#22713d')
axes[2,2].text(.04,.28,'Same four displacement samples, noise bounds, prior,\nand tolerances. Only their times change.',
               transform=axes[2,2].transAxes,fontsize=10,va='top',
               bbox={'facecolor':'#f5f6f7','edgecolor':'0.8','boxstyle':'round,pad=.5'})
fig.suptitle('Proper sampling restores hidden-state recovery in the discrete optimality system',fontsize=15)
fig.tight_layout(rect=(0,0,1,.97));fig.savefig(poor_root/'sampling_comparison.png',dpi=180);plt.close(fig)
print('Plotted the matched poor-versus-proper sampling comparison.')
