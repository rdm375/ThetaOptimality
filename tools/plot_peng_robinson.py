#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Render C++ outputs only: no EOS, inference, or qualification decisions."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
root=Path(__file__).resolve().parents[1]/'results'/'peng_robinson'
def rows(name):
    with (root/name).open() as f:return list(csv.DictReader(f))
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
q=rows('qualification.csv'); cells=rows('outer_cells.csv')
decisions={r['case']:r for r in rows('qualification_decisions.csv')}
cases=['one_density','four_densities','four_precise_densities']
titles=['One density\n0.1% error bound','Four densities\n0.1% error bounds','Four densities\n0.02% error bounds']
fig,axes=plt.subplots(1,3,figsize=(15,4.8),sharex=True,sharey=True)
for ax,case,title in zip(axes,cases,titles):
    selected=[r for r in cells if r['case']==case]
    polys=[]
    for r in selected:
        a,b,c,d=[100*float(r[k]) for k in ['ethane_lower','ethane_upper','nitrogen_lower','nitrogen_upper']]
        polys.append([(a,c),(b,c),(b,d),(a,d)])
    ax.add_collection(PolyCollection(polys,facecolors='#76a5d5',edgecolors='none',rasterized=True))
    ax.plot(10,5,'k*',ms=12,label='Synthetic truth')
    eq=[r for r in q if r['case']==case]; estimates={r['goal']:float(r['estimate']) for r in eq}
    ax.plot(100*estimates['ethane'],100*estimates['nitrogen'],'o',color='#b34128',label='Fitted representative')
    joint=decisions[case]['joint_qualified']=='1' # display the C++ joint verdict
    ax.set_title(title+'\n'+('QUALIFIED' if joint else 'NOT QUALIFIED'),color='#17643a' if joint else '#a33126')
    ax.set_xlim(0,16);ax.set_ylim(0,20);ax.set_xlabel('Ethane mole fraction (%)');ax.grid(alpha=.2)
axes[0].set_ylabel('Nitrogen mole fraction (%)');axes[0].legend(loc='upper right',fontsize=9)
fig.suptitle('Nonlinear outer sets from bounded density errors | all three tolerances: 1 percentage point',fontsize=13)
fig.text(.5,.015,'Blue: retained interval cells, including unresolved cells. One-density representative fixes nitrogen arbitrarily at 8%.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.04,1,.93));fig.savefig(root/'composition_sets.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(15,4.7))
for ax,case,title in zip(axes,cases,titles):
    records=[r for r in q if r['case']==case]
    for j,r in enumerate(records):
        lo,hi,est=[100*float(r[k]) for k in ['lower','upper','estimate']]
        color='#17643a' if r['qualified']=='1' else '#b34128'
        ax.plot([lo-est,hi-est],[j,j],lw=6,color=color,alpha=.65)
        ax.plot(0,j,'ko');ax.plot([-1,1],[j,j],'k|',ms=17)
    ax.set_yticks(range(3),[r['goal'].capitalize() for r in records]);ax.set_xlabel('Deviation from estimate (percentage points)');ax.set_title(title);ax.grid(axis='x',alpha=.2);ax.axvspan(-1,1,color='#cee6d4',alpha=.5,zorder=0)
fig.suptitle('C++ qualification verdict by component | bars: outer deviations; black ticks: ±1 percentage point tolerance',fontsize=13)
fig.tight_layout(rect=(0,0,1,.92));fig.savefig(root/'qualification.png',dpi=170);plt.close(fig)
obs=rows('observations.csv'); trace=rows('optimizer_trace.csv');surface=rows('objective_surface.csv')
fig,axes=plt.subplots(1,3,figsize=(15,4.7))
axes[0].errorbar(np.arange(1,5),[(float(r['fitted_density_kg_m3'])-float(r['observed_density_kg_m3']))/float(r['sigma_kg_m3']) for r in obs],yerr=1,fmt='o',capsize=5)
axes[0].axhline(0,color='k',lw=.8);axes[0].set_xticks(range(1,5));axes[0].set_xlabel('Observation index');axes[0].set_ylabel('(Fit − observation) / scale');axes[0].set_title('Four-density residuals\n0.1% scales')
E=np.array([float(r['ethane']) for r in surface]).reshape(101,101)*100
N=np.array([float(r['nitrogen']) for r in surface]).reshape(101,101)*100
C=np.array([float(r['cost']) for r in surface]).reshape(101,101)
axes[1].contour(E,N,C,levels=[1,2,5,20,100,1000],colors='#346e9e');axes[1].plot(10,5,'k*',ms=11);axes[1].set_xlabel('Ethane mole fraction (%)');axes[1].set_ylabel('Nitrogen mole fraction (%)');axes[1].set_title('Weighted least-squares objective\nLong valley = coupled uncertainties')
axes[2].semilogy([int(r['iteration']) for r in trace],[max(float(r['gradient_norm']),1e-15) for r in trace],'o-');axes[2].set_xlabel('Newton iteration');axes[2].set_ylabel('Adjoint gradient infinity norm');axes[2].set_title('Optimization convergence')
for ax in axes:ax.grid(alpha=.2)
fig.tight_layout();fig.savefig(root/'reconstruction.png',dpi=170);plt.close(fig)
print('Plotted C++ Peng–Robinson outputs; Python performs no EOS or inverse calculations.')
