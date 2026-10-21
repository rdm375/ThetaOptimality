# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Richard Myers

"""Plot C++ results. Python is an optional visualization dependency only."""
from pathlib import Path
import csv,json,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]
def read(name):
    with (root/'results'/name).open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
rows=read('convergence.csv');keys=['state','cost','tangent','adjoint','gradient','hessian','initial_state']
orders={}
fig,axes=plt.subplots(2,4,figsize=(15,7),layout='constrained')
for ax,key in zip(axes.flat,keys):
    for th in [.5,.65]:
        rr=[r for r in rows if abs(r['theta']-th)<1e-12]
        p=math.log2(rr[-2][key]/rr[-1][key]);orders.setdefault(str(th),{})[key]=p
        ax.loglog([r['h'] for r in rr],[r[key] for r in rr],'o-',label=f'θ={th:g}, p={p:.4f}')
    rr=[r for r in rows if r['theta']==.5];ax.loglog([r['h'] for r in rr],[rr[-1][key]*(r['h']/rr[-1]['h'])**2 for r in rr],'k--',alpha=.4,label='Slope 2')
    ax.set(title=key.replace('_',' ').capitalize(),xlabel='Maximum step h (min)',ylabel='Absolute error');ax.legend(fontsize=8);ax.grid(which='both',alpha=.2)
axes.flat[-1].axis('off');axes.flat[-1].text(.05,.8,'C++20 library validation\n\nSame forward theta scheme\nSame cost quadrature\nExact discrete tangent\nExact discrete adjoint\nIncremental-adjoint Hessian\n\nθ = 0.5: second order',va='top',fontsize=12)
fig.savefig(root/'results/convergence.png',dpi=150);plt.close(fig)
tr=read('trajectory.csv');fig,axes=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
for ax,computed,exact,title in zip(axes,['x','tlm_x','lambda_x'],['exact_x','exact_tlm_x','exact_lambda_x'],['Concentration x','Tangent dx/dξ','Cost-to-go adjoint λx']):
    ax.plot([r['time'] for r in tr],[r[exact] for r in tr],label='Analytical')
    ax.plot([r['time'] for r in tr[::8]],[r[computed] for r in tr[::8]],'o',ms=3,label='C++ theta, θ=0.5')
    ax.set(title=title,xlabel='Time t (min)',ylabel=title);ax.legend(fontsize=8);ax.grid(alpha=.2)
fig.savefig(root/'results/trajectory.png',dpi=150);plt.close(fig)
nt=read('newton.csv');fig,ax=plt.subplots(figsize=(6,4),layout='constrained')
ax.semilogy([r['iteration'] for r in nt],[max(r['error'],1e-16) for r in nt],'o-')
ax.set(title='Newton convergence at N=160, θ=0.5',xlabel='Iteration k',ylabel='Error to discrete minimizer |ξk − ξ*h|');ax.grid(alpha=.2)
fig.savefig(root/'results/newton.png',dpi=150);plt.close(fig)
controls=read('control_convergence.csv')
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for theta in [.5,.65]:
    rr=[r for r in controls if abs(r['theta']-theta)<1e-12]
    axes[0].loglog([r['h'] for r in rr[1:]],[r['absolute_error'] for r in rr[1:]],'o-',label=f'θ={theta:g}')
    axes[1].semilogx([r['h'] for r in rr],[r['discrete_minimizer'] for r in rr],'o-',label=f'θ={theta:g}')
axes[0].set(title='Convergence of optimized initial state',xlabel='Maximum step h (min)',ylabel='|ξ*h − ξ*|')
axes[0].grid(which='both',alpha=.2);axes[0].legend()
axes[1].axhline(controls[0]['analytical_minimizer'],color='k',ls='--',label='Analytical minimizer')
axes[1].set(title='Discrete control versus mesh size',xlabel='Maximum step h (min)',ylabel='Discrete minimizer ξ*h')
axes[1].grid(alpha=.2);axes[1].legend()
fig.savefig(root/'results/control_convergence.png',dpi=150);plt.close(fig)
(root/'results/orders.json').write_text(json.dumps(orders,indent=2))
table='| Quantity | θ = 0.5 | θ = 0.65 |\n|---|---:|---:|\n'
for k in keys:table+=f'| {k.replace("_"," ")} | {orders["0.5"][k]:.6f} | {orders["0.65"][k]:.6f} |\n'
(root/'results/orders.md').write_text(table)
print(table)
