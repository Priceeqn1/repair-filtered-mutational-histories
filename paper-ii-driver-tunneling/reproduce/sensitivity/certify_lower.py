"""Rational-residual parity witness for U=1/5, theta_L=1/1000, M=10.
Uses the original PDF model's source equilibrium, not frozen founder weights.
"""
from fractions import Fraction as F
from math import sqrt
import sys
# certification_base.py is included alongside this script
import certification_base as b

U=F(1,10)
b.TWO_L=(F(2),F(2),F(1,100))
b.bs=[b.ONE_C[0],b.ONE_L[0],b.TWO_C[0],b.TWO_L[0]]
b.ds=[b.ONE_C[1],b.ONE_L[1],b.TWO_C[1],b.TWO_L[1]]
b.cs=[F(0),F(0),b.TWO_C[2],b.TWO_L[2]]
b.THETA_MIN=F(1,100)

def source_bounds(nu):
    # x is positive root of (3-2U+4nu)x^2-(3-4U)x-4nu=0
    c0=3-2*U+4*nu;c1=3-4*U
    if nu==0:xl=xh=c1/c0
    else:
        xf=(float(c1)+sqrt(float(c1*c1+16*nu*c0)))/(2*float(c0))
        err=F(1,10**10)
        xl=F.from_float(xf)-err; xh=F.from_float(xf)+err
        p=lambda x:c0*x*x-c1*x-4*nu
        assert 0<xl<xh<1 and p(xl)<0<p(xh)
    ph=lambda x:(1+(2*U-1)*x)/(1+x)
    return ph(xh),ph(xl)

def certificate(M,nu):
    pre,post,Ip,Ik,Q,K,B,t=b.prepare(M,nu)
    qtilde=b.numeric_solve(K,t)
    qrat=[F.from_float(float(v)) for v in qtilde]
    rq=b.exact_residual(K,t,qrat)
    Eq=rq/b.THETA_MIN
    Bq=[sum((coeff*qrat[j] for j,coeff in row.items()),F(0)) for row in B]
    atilde=b.numeric_solve(Q,Bq)
    arat=[F.from_float(float(v)) for v in atilde]
    ra=b.exact_residual(Q,Bq,arat)
    Bnorm=max((sum(row.values(),F(0)) for row in B),default=F(0))
    Ea=F(M-1)/b.a*(ra+Bnorm*Eq)
    pl,ph=source_bounds(nu)
    ac,al=arat[Ip[(1,0)]],arat[Ip[(0,1)]]
    qc,ql=qrat[Ik[(0,0,1,0)]],qrat[Ik[(0,0,0,1)]]
    nlo=2*((1-ph)*(ac-Ea)+pl*(al-Ea))
    nhi=2*((1-pl)*(ac+Ea)+ph*(al+Ea))
    dlo=(1-ph)*(qc-Eq)+pl*(ql-Eq)
    dhi=(1-pl)*(qc+Eq)+ph*(ql+Eq)
    assert nlo>0 and dlo>0
    low,high=nlo/dhi,nhi/dlo
    print(f'M={M}, nu={nu} source phi [{float(pl):.12f},{float(ph):.12f}], numerical D interval [{float(low):.12f},{float(high):.12f}], bounds Eq={float(Eq):.3g}, Ea={float(Ea):.3g}')
    return low,high

if __name__=='__main__':
    results=[certificate(10,F(0)),certificate(10,F(1,5)),certificate(10,F(10))]
    assert results[0][0]>1 and results[1][1]<1 and results[2][0]>1
    print('CERTIFIED at least two finite M=10 parity crossings')
