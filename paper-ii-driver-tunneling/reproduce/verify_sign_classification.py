"""Reproduce exact unbounded sign identity and independent finite-M generator checks.
Requires: numpy, scipy.  Python 3.10+.  No external data or older scripts.
"""
import math
from itertools import product
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq

ONE_C=(306/755, 6/5)
ONE_L=(1/5, 3/2)
TWO_C=(12/5, 8/5, 1)
TWO_L=(2,2,1e-6)

def x_source(nu):
    # Positive equilibrium branch: (2+4nu)x^2 - x - 4nu=0.
    return (1+math.sqrt(1+16*nu*(2+4*nu)))/(2*(2+4*nu))

def phi_and_derivative(nu):
    x=x_source(nu)
    xp=4*(1-x*x)/(2*(2+4*nu)*x-1)
    return 1/(1+x), -xp/(1+x)**2

def success_q(b,d,theta):
    return 2*theta/(d+theta-b+math.sqrt((d+theta-b)**2+4*b*theta))

def q_damage(nu):
    B,d,theta=TWO_L
    Q=success_q(*TWO_C)
    k=d+theta-B+nu
    q=2*(theta+nu*Q)/(k+math.sqrt(k*k+4*B*(theta+nu*Q)))
    qp=(Q-q)/(2*B*q+k)
    return Q,q,qp

def unbounded(nu):
    bC,dC=ONE_C; b,dL=ONE_L
    f=bC/(dC-bC); a=dL-b; T=a*f-b
    Q,q,qp=q_damage(nu)
    phi,phip=phi_and_derivative(nu)
    h=a+nu; K=a*f*Q-b*q
    Z=(1-phi)*Q+phi*q
    W=f*nu*(Q-q)-T*q
    S=f*Q-phi*K/h
    D=2*S/Z
    derivative=2*(phi*K/(h*h*Z)+phip*Q*W/(h*Z*Z)-phi*Q*((1-phi)*T+f*nu)*qp/(h*Z*Z))
    return D,derivative,W

def compositions(n,k):
    if k==1:
        yield (n,); return
    for j in range(n+1):
        for tail in compositions(n-j,k-1):
            yield (j,)+tail

def finite_ratio(M,nu):
    if M<3: raise ValueError('Need M>=3 for sequential path')
    pre=[z for n in range(1,M) for z in compositions(n,2)]
    post=[z for n in range(1,M) for z in compositions(n,4)]
    ipre={z:i for i,z in enumerate(pre)}
    ipost={z:i for i,z in enumerate(post)}
    Q=lil_matrix((len(pre),len(pre)))
    K=lil_matrix((len(post),len(post)))
    B=lil_matrix((len(pre),len(post)))
    births=[ONE_C[0],ONE_L[0],TWO_C[0],TWO_L[0]]
    deaths=[ONE_C[1],ONE_L[1],TWO_C[1],TWO_L[1]]
    commits=[0,0,TWO_C[2],TWO_L[2]]
    def fill(states,index,G):
        t=np.zeros(len(states))
        for j,z in enumerate(states):
            exit_rate=0.
            for i,num in enumerate(z):
                if num==0: continue
                for intensity,delta in [(num*births[i],1),(num*deaths[i],-1)]:
                    exit_rate+=intensity
                    nz=list(z); nz[i]+=delta; nz=tuple(nz)
                    if 0<sum(nz)<M:
                        G[j,index[nz]]+=intensity
                if i in (1,3):
                    intensity=num*nu; exit_rate+=intensity
                    nz=list(z);nz[i]-=1;nz[i-1]+=1
                    G[j,index[tuple(nz)]]+=intensity
                if i>=2:
                    intensity=num*commits[i]; exit_rate+=intensity
                    t[j]+=intensity
            G[j,j]-=exit_rate
        return t
    fill(pre,ipre,Q)
    t=fill(post,ipost,K)
    for j,z in enumerate(pre):
        if sum(z)>=M-1: continue
        for i,num in enumerate(z):
            if num==0: continue
            # A mutation-marked birth adds a D daughter, preserving parent.
            destination=(z[0],z[1],int(i==0),int(i==1))
            B[j,ipost[destination]]+=num*births[i]
    q=spsolve(-K.tocsr(),t)
    A=spsolve(-Q.tocsr(),B.tocsr().dot(q))
    phi,_=phi_and_derivative(nu)
    s=(1-phi)*A[ipre[(1,0)]]+phi*A[ipre[(0,1)]]
    z=(1-phi)*q[ipost[(0,0,1,0)]]+phi*q[ipost[(0,0,0,1)]]
    return 2*s/z

def roots(fun,lo,hi):
    grid=np.r_[0,np.geomspace(lo,hi,180)]
    result=[]
    for x,y in zip(grid[:-1],grid[1:]):
        fx,fy=fun(x),fun(y)
        if fx*fy<0: result.append(brentq(fun,x,y,xtol=1e-12))
    return result

def main():
    for nu in [0,.001,.03,.25,2,10]:
        D,dD,W=unbounded(nu)
        step=1e-6 if 0<nu<.01 else (1e-7 if nu==0 else 1e-5)
        if nu==0:
            fd=(unbounded(step)[0]-D)/step
            assert abs(fd-dD)/(1+abs(dD))<.02
        else:
            fd=(unbounded(nu+step)[0]-unbounded(nu-step)[0])/(2*step)
            assert abs(fd-dD)/(1+abs(dD))<1e-6
        print(f'unbounded nu={nu:8g}: D={D:.9f}, Dprime={dD:+.9f}, W={W:+.6f}')
    print('Unbounded parity roots:',roots(lambda x:unbounded(x)[0]-1,1e-7,1e4))
    print('Unbounded stationary roots:',roots(lambda x:unbounded(x)[1],1e-7,1e4))
    print('Unbounded high-repair-level root:',roots(lambda x:unbounded(x)[2],1e-7,1e4))
    Q,_,_=q_damage(0)
    f=ONE_C[0]/(ONE_C[1]-ONE_C[0]);a=ONE_L[1]-ONE_L[0]
    T=a*f-ONE_L[0]
    print('Hazard certificate: f*B_L*Q-T =',f*TWO_L[0]*Q-T)
    assert f*TWO_L[0]*Q-T>0
    for M in [8,10,12]:
        data=[finite_ratio(M,v) for v in [0,.1,10,100]]
        print('finite M',M,'at 0,.1,10,100:',*[round(v,9) for v in data])
    assert finite_ratio(10,0)>1 and finite_ratio(10,.1)<1 and finite_ratio(10,10)>1
    assert finite_ratio(10,10)>finite_ratio(10,100)
    print('M=10 strict inequalities certify two parity crossings and two interior extrema.')

if __name__=='__main__':main()
