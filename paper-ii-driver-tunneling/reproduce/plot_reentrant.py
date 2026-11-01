"""Four-type reentrant pathway ratios. Run from the integrated project folder."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_sign_classification import unbounded,finite_ratio
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
p=Path(__file__).resolve().parents[1]/'figures'
p.mkdir(parents=True,exist_ok=True)
nu_un=np.concatenate(([0.0],np.geomspace(1e-7,300,200)))
xx=np.log10(1+nu_un)
fig,ax=plt.subplots(figsize=(6.75,3.65))
ax.plot(xx,[unbounded(float(t))[0] for t in nu_un],lw=2.1,label=r'$M=\infty$ (exact)')
nu_f=np.r_[0.0,np.geomspace(.0001,300,25)]
for M in (9,10):
    ys=[finite_ratio(M,float(t)) for t in nu_f]
    ax.plot(np.log10(1+nu_f),ys,linewidth=1.4,marker='o',markersize=2.3,label=f'$M={M}$ (finite generator)')
ax.axhline(1,linewidth=1,linestyle='--',color='gray')
ax.set_xlim(0,np.log10(301))
ax.set_xlabel(r'Log repair intensity $\log_{10}(1+\nu)$')
ax.set_ylabel(r'Sequential / direct ratio $D_M$')
ax.set_ylim(.75,1.115)
ax.grid(alpha=.15)
ax.legend(fontsize=8,frameon=False,ncol=2)
fig.tight_layout()
fig.savefig(p/'reentrant_pathways.pdf',bbox_inches='tight')
fig.savefig(p/'reentrant_pathways.png',dpi=150,bbox_inches='tight')
print('Saved figure; M9 and M10 values from finite generators.')
