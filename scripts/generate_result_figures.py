#!/usr/bin/env python3
"""Produce Figs. 2 and 3 of Paper I from their model equations.

Recreates the scientific data displayed in the manuscript's earlier PDFs;
line and marker styling may differ slightly from older figure snapshots.
"""
import csv
import sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from repair_histories.model import (two_hit_stats, event_distribution,
    survival_factor, equal_leading_burden_mutation_ratio)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=root/'figures'
DATA=root/'data'
OUT.mkdir(exist_ok=True)
DATA.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans',
                     'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})

def save_figure(fig, name):
    for ext in ['pdf','png']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

def two_hit_curve():
    # Manuscript: r=0.8, g=1.3, eps=0.5, nu in [0,2].
    nu_grid=np.linspace(0,2,201)
    d0,b0=two_hit_stats(0)
    direct=[]; abundance=[]
    for nu in nu_grid:
        d,b=two_hit_stats(float(nu))
        direct.append(d)
        abundance.append(b/b0)
    fig,ax=plt.subplots(figsize=(7.3,4.3))
    ax.plot(nu_grid,abundance,label='Normalized two-hit abundance',lw=2.4)
    ax.plot(nu_grid,direct,label='Direct-route share',lw=2.4,ls='--')
    ax.set(xlabel=r'Repair rate $\nu$',ylabel='Leading-order fraction',xlim=(0,2),ylim=(0,1.02))
    ax.grid(alpha=.18)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig,'paper_I_two_hit_mechanism')
    with (DATA/'figure_2_two_hit_curve.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['nu','normalized_two_hit_abundance','direct_route_share'])
        w.writerows(zip(nu_grid,abundance,direct))
    assert abs(abundance[100]-.08118811881188119) < 1e-5
    assert abs(direct[100]-25/41) < 1e-12
    print('Figure 2: nu=1 burden ratio %.8f, direct share %.8f'%(abundance[100],direct[100]))

def equal_burden_event_counts():
    # Manuscript: fixed k=7; leading burden matching U2/U1 = [E7(0)/E7(1)]^(1/7).
    k=7
    p0=event_distribution(k,0);p1=event_distribution(k,1)
    e0=survival_factor(k,0);e1=survival_factor(k,1)
    ratio=equal_leading_burden_mutation_ratio(k,0,1)
    x=np.arange(1,k+1)
    fig,ax=plt.subplots(figsize=(7.3,4.4))
    ax.plot(x,p0,'o-',lw=2.2,label=rf'$\nu=0$ ($E_7={e0:.3g}$)')
    ax.plot(x,p1,'s-',lw=2.2,label=rf'$\nu=1$ ($E_7={e1:.4g}$)')
    ax.set(xlabel=r'Number of mutational jumps $J$',
           ylabel=r'Conditional auxiliary-route probability $\Pr(J\mid k=7)$',
           xticks=x, xlim=(.7,k+.3),ylim=(0,.46))
    ax.grid(axis='y',alpha=.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig,'paper_I_jump_count')
    with (DATA/'figure_3_equal_burden_events.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['jump_count','prob_nu0','prob_nu1'])
        w.writerows(zip(x,p0,p1))
    assert abs(e0-1)<1e-10
    assert abs(e1-.00956703934)<2e-8
    assert abs(ratio-1.94294429)<2e-6
    print('Figure 3: E7(1)=%.11f, matched leading U2/U1=%.8f'%(e1,ratio))

if __name__=='__main__':
    two_hit_curve()
    equal_burden_event_counts()
