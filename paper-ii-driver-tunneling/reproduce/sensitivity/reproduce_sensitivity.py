"""Sensitivity analysis of exact unbounded and finite-M repair-dependent pathway ratios.
Python 3.10+, numpy scipy matplotlib.
All time rates are dimensionless with pristine resident division as time unit.
No clinical calibration is assumed.
"""
import math
from itertools import product
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

SOURCE_R=0.5
SOURCE_G=2.0
SOURCE_EPS=0.25
ONE_C=(306/755, 6/5)
ONE_L=(1/5,3/2)
TWO_C=(12/5,8/5,1.)
TWO_L_BIRTH=2.

def resident_x(U,nu):
    """Physical equilibrium for U<3/4: (3-2U+4nu)x²-(3-4U)x-4nu=0."""
    U=np.asarray(U);nu=np.asarray(nu)
    a=3-2*U+4*nu
    b=3-4*U
    return (b+np.sqrt(b*b+16*nu*a))/(2*a)

def founder_phi(U,nu):
    x=resident_x(U,nu)
    return (1+(2*U-1)*x)/(1+x)

def success_q(b,d,theta):
    # stable expression for positive solution of b q²+(d+theta-b)q-theta=0
    return 2*theta/(d+theta-b+np.sqrt((d+theta-b)**2+4*b*theta))

def unbounded(nu,U=.5,theta_L=1e-6,delta_L=2.,b_L_2=2.,f_C=None):
    nu=np.asarray(nu)
    if f_C is None:
        f_C=ONE_C[0]/(ONE_C[1]-ONE_C[0])
    Q=success_q(*TWO_C)
    k=delta_L+theta_L-b_L_2+nu
    q=2*(theta_L+nu*Q)/(k+np.sqrt(k*k+4*b_L_2*(theta_L+nu*Q)))
    bL,dL=ONE_L
    AL=(bL*q+nu*f_C*Q)/(dL-bL+nu)
    phi=founder_phi(U,nu)
    D=2*((1-phi)*f_C*Q+phi*AL)/((1-phi)*Q+phi*q)
    return D

def finite_ratio(M,nu,U=.5,theta_L=1e-6,delta_L=2.,b_L_2=2.,one_c=ONE_C):
    if M<3:raise ValueError('M>=3 required')
    pre=[z for n in range(1,M) for z in compositions(n,2)]
    post=[z for n in range(1,M) for z in compositions(n,4) if z[2]+z[3]>0]
    ip={z:i for i,z in enumerate(pre)}; ik={z:i for i,z in enumerate(post)}
    Q=lil_matrix((len(pre),len(pre)));K=lil_matrix((len(post),len(post)));B=lil_matrix((len(pre),len(post)))
    br=[one_c[0],ONE_L[0],TWO_C[0],b_L_2]
    dr=[one_c[1],ONE_L[1],TWO_C[1],delta_L]
    tr=[0,0,TWO_C[2],theta_L]
    def fill(states,indices,MAT):
        t=np.zeros(len(states))
        for i,z in enumerate(states):
            total=sum(z);out=0.
            for k,n in enumerate(z):
                if not n: continue
                for r,change in [(n*br[k],1),(n*dr[k],-1)]:
                    out+=r
                    nz=list(z);nz[k]+=change;nz=tuple(nz)
                    j=indices.get(nz)
                    if j is not None: MAT[i,j]+=r
                if k in (1,3):
                    r=n*nu;out+=r
                    nz=list(z);nz[k]-=1;nz[k-1]+=1
                    MAT[i,indices[tuple(nz)]]+=r
                if k>=2:
                    r=n*tr[k];out+=r;t[i]+=r
            MAT[i,i]-=out
        return t
    fill(pre,ip,Q)
    t=fill(post,ik,K)
    for i,z in enumerate(pre):
        if sum(z)>=M-1:continue
        for k,n in enumerate(z):
            if not n:continue
            nz=list(z)+[0,0];nz[2+k]+=1
            B[i,ik[tuple(nz)]]+=n*br[k]
    q=spsolve(-K.tocsr(),t)
    A=spsolve(-Q.tocsr(),B.tocsr().dot(q))
    phi=float(founder_phi(U,nu))
    return 2*((1-phi)*A[ip[(1,0)]]+phi*A[ip[(0,1)]])/((1-phi)*q[ik[(0,0,1,0)]]+phi*q[ik[(0,0,0,1)]])

def compositions(n,k):
    if k==1:
        yield (n,);return
    for i in range(n+1):
        for tail in compositions(n-i,k-1):yield (i,)+tail

def min_unbounded(U,theta_L=1e-6,f_C=None):
    obj=lambda y:float(unbounded(10**y,U=U,theta_L=theta_L,f_C=f_C))
    s=minimize_scalar(obj,bounds=(-9,6),method='bounded',options={'xatol':1e-9})
    vals=[(float(unbounded(0,U=U,theta_L=theta_L,f_C=f_C)),0),(float(s.fun),10**s.x),(obj(6),1e6)]
    return min(vals)

def unbounded_parity_roots(U,theta_L=1e-6,f_C=None):
    ns=np.r_[0,np.geomspace(1e-9,1e6,600)]
    y=unbounded(ns,U=U,theta_L=theta_L,f_C=f_C)-1
    return [brentq(lambda v:float(unbounded(v,U=U,theta_L=theta_L,f_C=f_C)-1),a,b,xtol=1e-11)
            for a,b,u,v in zip(ns[:-1],ns[1:],y[:-1],y[1:]) if u*v<0]

if __name__=='__main__':
    # source validation; original manuscript U=.5 must give phi=.66667 at zero repair
    assert abs(float(founder_phi(.5,0))-2/3)<1e-12
    assert abs(float(founder_phi(.2,0))-4/15)<1e-12
    # verify source formula for old example at U=.5
    assert abs(float(unbounded(0,U=.5))-1.0183398951616112)<1e-10
    for U,theta in [(.5,1e-6),(.2,1e-3),(.1,1e-2)]:
        a,b=min_unbounded(U,theta)
        print('Unbounded:',U,theta,'D0',float(unbounded(0,U,theta)),'minimum',a,'nu',b,'roots',unbounded_parity_roots(U,theta))
        print('Finite M=10:',*[round(finite_ratio(10,r,U,theta),9) for r in [0,.1,.2,10]])
