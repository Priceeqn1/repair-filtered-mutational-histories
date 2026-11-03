"""Journal figure for the Paper II finite-M=10 robustness comparison.
Curves: numerical inversion of the full killed generators (not certified across
all repair values). Colored markers indicate repair rates at which strict parity
inequalities were checked separately with rational-residual certificates.
No empirical calibration is implied.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

SENS=Path(__file__).resolve().parent/'sensitivity'
sys.path.insert(0,str(SENS))
from reproduce_sensitivity import finite_ratio, unbounded
OUT=Path(__file__).resolve().parent.parent/'figures';OUT.mkdir(exist_ok=True)

cases=[(.5,1e-6,'Original existence example','tab:blue',0.1),(.2,1e-3,'Moderate source/commitment','tab:orange',0.1),(.1,1e-2,'Lower lesion input','tab:green',0.2)]
nu=np.geomspace(1e-4,1e2,145)
TH={1e-6:r'10^{-6}',1e-3:'0.001',1e-2:'0.01'}
fig,ax=plt.subplots(figsize=(9.0,5.3),constrained_layout=True)
for U,theta,label,color,middle in cases:
    D=np.array([finite_ratio(10,float(v),U=U,theta_L=theta) for v in nu])
    ax.plot(nu,D,lw=2.1,color=color,label=rf'$U={U:g}$, $\theta_L={TH[theta]}$')
    witnesses=[middle,10]
    y=[finite_ratio(10,float(v),U=U,theta_L=theta) for v in witnesses]
    ax.scatter(witnesses,y,s=39,color=color,edgecolor='white',zorder=5,lw=.6)
ax.axhline(1,color='0.35',lw=1.2,ls='--')
ax.set_xscale('log');ax.set_xlim(1e-4,1e2);ax.set_ylim(.946,1.034)
ax.set_xlabel(r'Repair intensity $\nu$ (log scale)',fontsize=12)
ax.set_ylabel(r'Sequential / direct ratio $D_{10}$',fontsize=12)
ax.legend(loc='lower left',fontsize=10,frameon=False)
ax.text(1.8e-4,1.031,r'$D_{10}>1$: sequential favored',fontsize=9,color='0.34',va='top')
ax.text(1.8e-4,0.988,r'$D_{10}<1$: direct favored',fontsize=9,color='0.34')
ax.tick_params(labelsize=10)
ax.minorticks_on()
ax.spines[['top','right']].set_visible(False)
# Focus inset on lower-U nearly-parity effect, with independent evaluated points.
axins=ax.inset_axes([.70,.14,.28,.27])
local=np.linspace(.045,.6,80)
near=np.array([finite_ratio(10,float(v),U=.1,theta_L=.01) for v in local])
axins.plot(local,near,color='tab:green',lw=1.8)
axins.axhline(1,color='0.4',ls='--',lw=.9)
axins.scatter([.2],[finite_ratio(10,.2,U=.1,theta_L=.01)],color='tab:green',edgecolor='white',s=45,zorder=5)
axins.set_xlim(.045,.6);axins.set_ylim(.9991,1.0017)
axins.set_title(r'$U=0.1,\ \theta_L=0.01$: near parity',fontsize=9,pad=3)
axins.tick_params(axis='both',labelsize=8)
axins.set_xlabel(r'$\nu$ (linear)',fontsize=8)
axins.ticklabel_format(axis='y',useOffset=False)
for ext in ('pdf','png'):
    fig.savefig(OUT/f'finite_M10_sensitivity.{ext}',dpi=230,bbox_inches='tight')
plt.close(fig)
# supplementary source-only sweep, fast unbounded analytic expressions
fig2,ax2=plt.subplots(figsize=(8.0,4.5),constrained_layout=True)
nu2=np.geomspace(1e-5,1e2,250)
for U in [.5,.2,.1,.01]:
    ax2.plot(nu2,unbounded(nu2,U=U,theta_L=1e-6),label=rf'$U={U:g}$',lw=1.9)
ax2.axhline(1,ls='--',color='0.4',lw=1)
ax2.set_xscale('log');ax2.set_xlabel(r'Repair intensity $\nu$ (log scale)')
ax2.set_ylabel(r'Sequential / direct ratio $D_\infty$')
ax2.legend(frameon=False);ax2.spines[['top','right']].set_visible(False)
for ext in ('pdf','png'):
    fig2.savefig(OUT/f'lesion_input_only.{ext}',dpi=210,bbox_inches='tight')
plt.close(fig2)
print('Created finite_M10_sensitivity.pdf and lesion_input_only.pdf')
