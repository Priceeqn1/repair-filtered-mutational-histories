"""Figure: certified-example finite-M ratios, using same generator as checks."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_sign_classification import finite_ratio
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
Ms=list(range(3,12))
ns=[0.0,0.1,10.0,1000.0]
y={nu:[finite_ratio(m,nu) for m in Ms] for nu in ns}
fig,ax=plt.subplots(figsize=(6.4,3.65))
styles=['-o','--s',':^','-.d']
for nu,st in zip(ns,styles):
    ax.plot(Ms,y[nu],st,markersize=3.7,lw=1.6,label=rf'$\nu={nu:g}$')
ax.axhline(1,ls='--',lw=1.0,color='.35',label='equal pathways')
ax.set_xlabel(r'Establishment threshold $M$ (cells)')
ax.set_ylabel(r'Sequential / direct ratio $D_M$')
ax.set_xlim(2.8,11.3)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_ylim(.66,1.08)
ax.legend(frameon=False,ncol=2,fontsize=8)
ax.grid(alpha=.14)
fig.tight_layout()
fig.savefig(Path(__file__).resolve().parents[1]/'figures'/'threshold_regimes.pdf',bbox_inches='tight')
fig.savefig(Path(__file__).resolve().parents[1]/'figures'/'threshold_regimes.png',dpi=160,bbox_inches='tight')
print('M  D(0) D(.1) D(10) D(1000)')
for k,m in enumerate(Ms):print(m,*[round(y[n][k],7) for n in ns])
