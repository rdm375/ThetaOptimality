#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot retained C++ two-tank results; no simulation or qualification logic."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]/'results/two_tank'
def rows(name):
    with (root/name).open() as f:return list(csv.DictReader(f))
def col(r,k):return np.array([float(x[k]) for x in r])
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
r=rows('trajectories.csv');obs=rows('observations.csv');t=col(r,'elapsed_min')
fig,axes=plt.subplots(2,2,figsize=(12,8))
for ax,key,title in [(axes[0,0],'c1','Small downstream tank'),(axes[0,1],'c2','Large upstream tank')]:
    ax.plot(t,col(r,'true_'+key+'_mg_L'),color='0.65',lw=2.5,label='Synthetic truth')
    ax.plot(t,col(r,'analytical_fit_'+key+'_mg_L'),'k--',label='Analytical noisy-data fit')
    ax.plot(t,col(r,'discrete_'+key+'_mg_L'),color='#2876a5',label='Discrete fit, θ=0.5, N=160')
    ax.set_title(title);ax.set_xlabel('Elapsed time from t₀ (min)');ax.set_ylabel('Concentration (mg/L)');ax.grid(alpha=.2)
axes[0,0].errorbar(col(obs,'elapsed_min'),col(obs,'observed_c1_mg_L'),yerr=col(obs,'sigma_mg_L'),fmt='o',color='#ab4e30',capsize=3,label='Downstream data ± scale');axes[0,0].legend(fontsize=8)
axes[1,0].plot(t,col(r,'true_inventory_mg'),color='0.65',label='Truth');axes[1,0].plot(t,col(r,'discrete_inventory_mg'),color='#2876a5',label='Discrete fitted state');axes[1,0].set_ylabel('Total tracer inventory (mg)');axes[1,0].set_xlabel('Elapsed time (min)');axes[1,0].legend();axes[1,0].grid(alpha=.2)
for key in ['c1','c2']:axes[1,1].plot(t,1000*(col(r,'discrete_'+key+'_mg_L')-col(r,'analytical_fit_'+key+'_mg_L')),label=key+' difference')
axes[1,1].set_xlabel('Elapsed time (min)');axes[1,1].set_ylabel('Discrete fit − analytical fit (μg/L)');axes[1,1].legend();axes[1,1].grid(alpha=.2)
fig.suptitle('Two-tank state finding | known flow and volumes | downstream-only baseline observations')
fig.tight_layout(rect=(0,0,1,.95));fig.savefig(root/'reconstruction.png',dpi=170);plt.close(fig)
gain=rows('amplification.csv');worst=rows('worst_direction.csv');summary={x['quantity']:float(x['value']) for x in rows('summary.csv')}
fig,axes=plt.subplots(1,3,figsize=(15,4.7))
axes[0].plot(col(gain,'elapsed_min'),col(gain,'analytical_gain'),'k--',label='Analytical gain');axes[0].plot(col(gain,'elapsed_min'),col(gain,'discrete_gain'),color='#2876a5',label='Discrete tangent-map gain');axes[0].axhline(1,color='0.5',lw=.7);axes[0].axvline(summary['analytical_peak_time_min'],color='#ab4e30',ls=':',label=f"Peak {summary['analytical_peak_gain']:.6f} at {summary['analytical_peak_time_min']:.6f} min");axes[0].set_xlabel('Elapsed time (min)');axes[0].set_ylabel('Maximum Euclidean concentration gain');axes[0].set_title('Stable eigenvalues; transient amplification');axes[0].legend(fontsize=8)
axes[1].plot(col(worst,'elapsed_min'),col(worst,'perturbation_c1'),label='Downstream perturbation');axes[1].plot(col(worst,'elapsed_min'),col(worst,'perturbation_c2'),label='Upstream perturbation');axes[1].plot(col(worst,'elapsed_min'),col(worst,'euclidean_gain'),'k--',label='Euclidean norm');axes[1].plot(col(worst,'elapsed_min'),col(worst,'weighted_mass_error_bound')/float(worst[0]['weighted_mass_error_bound']),':',color='#287b43',label='Perturbation mass / initial mass');axes[1].set_xlabel('Elapsed time (min)');axes[1].set_ylabel('Normalized tangent perturbation / gain');axes[1].set_title('Worst initial direction, propagated by the TLM');axes[1].legend(fontsize=8)
axes[2].plot(col(worst,'perturbation_c1'),col(worst,'perturbation_c2'),color='#2876a5');axes[2].plot(float(worst[0]['perturbation_c1']),float(worst[0]['perturbation_c2']),'ko',mfc='white',label='Initial direction');axes[2].plot(float(worst[-1]['perturbation_c1']),float(worst[-1]['perturbation_c2']),'ks',label='End');axes[2].set_xlabel('Normalized downstream perturbation');axes[2].set_ylabel('Normalized upstream perturbation');axes[2].set_title('Perturbation phase plane');axes[2].legend(fontsize=8)
for ax in axes:ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(root/'amplification.png',dpi=170);plt.close(fig)
conv=rows('convergence.csv');orders={(float(x['theta']),x['component']):float(x['order']) for x in rows('orders.csv')};keys=['state_c1','state_c2','tangent','nodal_adjoint','control_c2','nodal_cost','nodal_gradient','nodal_hessian','running_cost','running_gradient','running_costate','running_hessian']
fig,axes=plt.subplots(3,4,figsize=(15,10))
for ax,key in zip(axes.flat,keys):
    for theta,marker in [(.5,'o'),(.65,'s')]:
        data=[x for x in conv if float(x['theta'])==theta];ax.loglog(col(data,'h_min'),col(data,key),marker+'-',ms=4,label=f'θ={theta}; order {orders[(theta,key)]:.3f}')
    ax.set_title(key.replace('_',' '));ax.set_xlabel('Step size (min)');ax.set_ylabel('Normalized error');ax.grid(which='both',alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Analytical temporal-convergence comparisons | directly observed c₁(t₀) control is roundoff-limited')
fig.tight_layout(rect=(0,0,1,.95));fig.savefig(root/'convergence.png',dpi=170);plt.close(fig)
qd=rows('qualification_decisions.csv');qe=rows('qualification_envelopes.csv');cases=['bounded_primary','bounded_refined','added_upstream','joint95_refined'];labels=['Bounded data | N=80','Bounded data | N=1280','Added upstream sensor | N=1280','Joint 95% Gaussian noise box | N=1280'];goals=['downstream_concentration','upstream_concentration','inventory'];titles=['Downstream | tolerance 0.005 mg/L','Upstream | tolerance 0.012 mg/L','Inventory | tolerance 0.04 mg']
fig,axes=plt.subplots(4,3,figsize=(15,13),sharex=True)
for i,case in enumerate(cases):
    for j,goal in enumerate(goals):
        ax=axes[i,j];data=[x for x in qe if x['case']==case and x['goal']==goal];d=next(x for x in qd if x['case']==case and x['goal']==goal);tol=float(d['tolerance']);whole=[x for x in data if x['whole_bound']!=''];bound=col(whole,'whole_bound')/tol
        ax.step(np.r_[col(whole,'elapsed_min'),float(whole[-1]['interval_end_min'])],np.r_[bound,bound[-1]],where='post',color='#b45b27',label='Physical whole-interval bound');ax.plot(col(data,'elapsed_min'),col(data,'analytical_bound')/tol,'--',color='#2876a5',label='Analytical physical-set reference');ax.plot(col(data,'elapsed_min'),col(data,'local95_indicator')/tol,':',color='0.6',label='Gaussian ellipse indicator only');ax.axhline(1,color='k',lw=.8)
        if d['first_below_min'] and d['first_below_min']!=d['persistent_min']:ax.axvline(float(d['first_below_min']),color='#b45b27',ls=':',label='First below, not persistent')
        status='NOT QUALIFIED'
        if d['qualified']=='1':status='QUALIFIED from '+d['persistent_min']+' min';ax.axvline(float(d['persistent_min']),color='#287b43',ls='--')
        ax.set_title(titles[j]+'\n'+status,fontsize=9,color='#287b43' if d['qualified']=='1' else '#a13c2a');ax.set_yscale('log');ax.set_ylim(.02,max(2.,1.4*np.max(bound)));ax.grid(which='both',alpha=.15)
        if j==0:ax.set_ylabel(labels[i]+'\nError bound / tolerance')
        if i==3:ax.set_xlabel('Elapsed time (min)')
axes[1,0].legend(fontsize=7);fig.suptitle('Persistent engineering qualification | theta = 0.5 | whole intervals through 10 min');fig.text(.5,.01,'All-data retrospective decisions. Deterministic bounded errors and Gaussian noise boxes are different assumptions.',ha='center');fig.tight_layout(rect=(0,.025,1,.96));fig.savefig(root/'qualification.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(13,4.8));ax.axis('off');cells=[];colors=[]
for case,label in zip(cases,labels):
    ds=[next(x for x in qd if x['case']==case and x['goal']==g) for g in goals];row=[label];color=['#f1f3f5']
    for d in ds:row.append('From '+d['persistent_min']+' min' if d['qualified']=='1' else 'NOT QUALIFIED');color.append('#d8ecdf' if d['qualified']=='1' else '#f4deda')
    d=ds[0];row.append('From '+d['joint_min']+' min' if d['joint_qualified']=='1' else 'NOT QUALIFIED');color.append('#d8ecdf' if d['joint_qualified']=='1' else '#f4deda');cells.append(row);colors.append(color)
table=ax.table(cellText=cells,colLabels=['Experiment','Downstream','Upstream','Inventory','ALL GOALS'],cellColours=colors,colWidths=[.32,.17,.17,.17,.17],cellLoc='center',loc='center');table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,2.6);ax.set_title('TWO-TANK ENGINEERING QUALIFICATION\n0.005 mg/L downstream | 0.012 mg/L upstream | 0.04 mg inventory',pad=18);fig.text(.5,.04,'Elapsed times; numerical error included; every remaining whole interval must satisfy the goals.',ha='center');fig.tight_layout(rect=(0,.08,1,1));fig.savefig(root/'qualification_verdict.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4.4))
for ax,key in zip(axes,['c1','c2']):ax.plot(t,col(r,'exact_probe_running_costate_'+key),'k--',label='Analytical');ax.plot(t,col(r,'discrete_probe_running_costate_'+key),color='#2876a5',label='Discrete cost-to-go convention');ax.set_xlabel('Elapsed time (min)');ax.set_ylabel('Costate ((mg/L) min)');ax.set_title('Running-cost '+key+' costate');ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(root/'costate.png',dpi=170);plt.close(fig)
print('Rendered C++ two-tank results; Python performs no model, optimization, or qualification calculations.')
