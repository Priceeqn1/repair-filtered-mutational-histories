"""Generate the exact five-state opportunity-establishment illustration."""
import itertools
from pathlib import Path
import numpy as np
from scipy.linalg import solve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
STATES=[z for n in (1,2) for z in itertools.product(range(n+1),repeat=2) if sum(z)==n]
I={z:i for i,z in enumerate(STATES)}

def vals(nu):
    Q=np.zeros((len(STATES),len(STATES)))
    b=np.zeros(len(STATES));e=np.zeros(len(STATES))
    for i,z in enumerate(STATES):
        for k,(birth,death) in enumerate([(10,1),(4,1)]):
            if not z[k]: continue
            rate=birth*z[k];Q[i,i]-=rate
            if sum(z)==2: e[i]+=rate
            else:
                zz=list(z);zz[k]+=1;Q[i,I[tuple(zz)]]+=rate;b[i]+=rate
            rate=death*z[k];Q[i,i]-=rate
            if sum(z)==2:
                zz=list(z);zz[k]-=1;Q[i,I[tuple(zz)]]+=rate
        Q[i,i]-=nu*z[1]
        if z[1]:
            zz=list(z);zz[1]-=1;zz[0]+=1;Q[i,I[tuple(zz)]]+=nu*z[1]
    return (solve(-Q,b)[I[(0,1)]],solve(-Q,e)[I[(0,1)]])

grid=np.r_[np.linspace(0,.5,120),np.linspace(.5,20,300)]
y=np.array([vals(z) for z in grid])
fig,ax=plt.subplots(figsize=(7,4))
l1=ax.plot(grid,y[:,0],lw=2,label='Eligible births')[0]
ax.set_xlabel('Repair rate $\\nu$'); ax.set_ylabel('Expected eligible births')
ax2=ax.twinx()
l2=ax2.plot(grid,y[:,1],linestyle='--',lw=1.7,label='Establishment probability')[0]
ax2.set_ylabel('Probability of establishing (size 3)')
ax.set_title('Mutation opportunity vs. premature establishment ($M=3$)')
ax.legend(handles=[l1,l2],loc='center right',frameon=False)
fig.tight_layout()
out = Path(__file__).resolve().parents[1] / 'figures'
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / 'opportunity_tradeoff_M3.pdf',bbox_inches='tight')
fig.savefig(out / 'opportunity_tradeoff_M3.png',dpi=180,bbox_inches='tight')
print('plot generated; opportunity minimum',grid[np.argmin(y[:,0])],min(y[:,0]))
