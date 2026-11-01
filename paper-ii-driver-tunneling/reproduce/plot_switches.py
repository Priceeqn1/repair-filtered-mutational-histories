"""Generate population-weighted and conditional three-cell M=3 pathway curves.

Model functions come from the independent finite-generator verification.
No arbitrary calibrated rates: all rates are the exact mathematical examples.
"""
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify_generator import model, pop

OUT=Path(__file__).resolve().parents[1]/'figures'
OUT.mkdir(parents=True,exist_ok=True)
nu=np.linspace(0,3.5,260)
conditional=np.array([model(float(v))[0] for v in nu])
weighted=np.array([pop(float(v)) for v in nu])
nu_pop=brentq(lambda x: pop(x)-1,0,1,xtol=1e-14)
nu_cond=brentq(lambda x: model(x)[0]-1,0,3.5,xtol=1e-14)
assert abs(nu_pop-0.192345334861229)<1e-10
assert abs(nu_cond-1.10964042855)<1e-9
fig,ax=plt.subplots(figsize=(7.0,4.25))
ax.plot(nu,weighted,label='Exact population-weighted ratio',lw=2.0)
ax.plot(nu,conditional,':',label='Conditional lesion-bearing ratio',lw=1.8)
ax.axhline(1,color='.45',ls='--',lw=1)
ax.axvline(nu_pop,color='.6',ls=':',lw=.8)
ax.scatter([nu_pop],[1],s=32,zorder=4)
ax.annotate(r'$\nu_c^{\mathrm{pop}}\approx0.19235$',xy=(nu_pop,1),xytext=(.35,1.53),arrowprops=dict(arrowstyle='-',lw=.8),fontsize=9)
ax.set_xlim(0,3.5)
ax.set_ylim(0,1.75)
ax.set_xlabel(r'Repair rate $\nu$')
ax.set_ylabel('Sequential / direct pathway ratio')
ax.legend(loc='lower right',fontsize=8,frameon=False)
fig.tight_layout()
fig.savefig(OUT/'population_repair_switch_M3.pdf',bbox_inches='tight')
fig.savefig(OUT/'population_repair_switch_M3.png',dpi=180,bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(figsize=(6.6,3.85))
ax.plot(nu,conditional,lw=2,label=r'$D_3(\nu)$, lesion-bearing founder')
ax.axhline(1,color='.4',ls='--',lw=1,label='Equal pathways')
ax.axvline(nu_cond,color='.65',ls=':',lw=.8)
ax.scatter([nu_cond],[1],s=30,zorder=4)
ax.set_xlabel(r'Repair rate $\nu$')
ax.set_ylabel('Sequential / direct pathway ratio')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT/'repair_switch_M3.pdf',bbox_inches='tight')
fig.savefig(OUT/'repair_switch_M3.png',dpi=180,bbox_inches='tight')
plt.close(fig)
print('Created population_repair_switch_M3.pdf and supplementary repair_switch_M3.pdf')
