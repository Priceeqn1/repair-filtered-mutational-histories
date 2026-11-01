"""Independent numeric and exact rational checks for the M=3 C,L,D model."""
import itertools
import numpy as np
from fractions import Fraction
from scipy.optimize import brentq
from scipy.linalg import solve

BIRTH=np.array([4.0,.1,.5]); DEATH=np.ones(3); THETA=50.0
STATES=[z for n in (1,2) for z in itertools.product(range(n+1),repeat=3) if sum(z)==n]
IDX={z:i for i,z in enumerate(STATES)}

def matrices(nu):
    n=len(STATES);Q=np.zeros((n,n));P=np.zeros((n,n));t=np.zeros(n)
    for i,z in enumerate(STATES):
        z=np.array(z); size=sum(z)
        t[i]=THETA*z[2]
        Q[i,i]-=t[i]
        for k in range(3):
            if z[k] == 0: continue
            rate=BIRTH[k]*z[k]; Q[i,i]-=rate
            if size+1<3:
                w=z.copy();w[k]+=1;Q[i,IDX[tuple(w)]]+=rate
                if k<2:
                    w2=z.copy();w2[2]+=1
                    P[i,IDX[tuple(w2)]]+=rate;P[i,IDX[tuple(w)]]-=rate
            rate=DEATH[k]*z[k];Q[i,i]-=rate
            if size-1>0:
                w=z.copy();w[k]-=1;Q[i,IDX[tuple(w)]]+=rate
        rate=nu*z[1]
        Q[i,i]-=rate
        if rate:
            w=z.copy();w[1]-=1;w[0]+=1;Q[i,IDX[tuple(w)]]+=rate
    return Q,P,t

def model(nu):
    Q,B,t=matrices(nu);p=solve(-Q,t);dp=solve(-Q,B@p)
    q=p[IDX[(0,0,1)]]
    ac=dp[IDX[(1,0,0)]];al=dp[IDX[(0,1,0)]]
    return 2*al/q,2*ac/q,q

def formula(nu):
    P=40000*nu**4+2397000*nu**3+15730820*nu**2+16064982*nu+1592283
    Q=1000*nu**3+8300*nu**2+14430*nu+6771
    return 2575*P/(11752*(5*nu+263)*Q)

def xof(nu): return (1+np.sqrt(64*nu**2+32*nu+1))/(8*nu+4)

def pop(nu):
    dl,dc,_=model(nu); x=xof(nu)
    return (x*dc+dl)/(1+x)

for nu in (0.,.1,.192345334861229,1.,1.10964042855,3.,8.):
    dl,dc,q=model(nu)
    assert np.isclose(dl,formula(nu),rtol=1e-11), (nu,dl,formula(nu))
    assert np.isclose(q,10400/10607,rtol=1e-11)
    assert np.isclose(dc,2575/1469,rtol=1e-11)
    print(f'nu={nu:.12g} DL={dl:.12g} DC={dc:.12g} population={pop(nu):.12g} q={q:.12g}')
root=brentq(lambda v:pop(v)-1,0,2,xtol=1e-14)
print('root',root)
assert abs(root-0.192345334861229)<1e-11
print('ALL MODEL-GENERATOR CHECKS PASS')
