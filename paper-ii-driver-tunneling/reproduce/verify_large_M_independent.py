"""Independent finite-state audit for the subcritical large-M theorem.

Post-hit state (C,L,D): at least one D alive; death of last D is killing.
At each ordinary one-driver birth, an inserted D daughter is counted at first
order in rare V; an insertion from size M-1 is excluded.
"""
from itertools import product
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

BC,BL,BD = 0.375, 0.095, 0.5
DC,DL,DD = 1., 2., 1.
THETA = 1.
U,R,G,EPS = .5,.5,2.,.25


def founder_fraction(nu):
    # Exact physical equilibrium of the one-locus resident for these parameters.
    x=(1+np.sqrt(64*nu*nu+32*nu+1))/(8*nu+4)
    return (U*x+R*(1-x))/(x+R*(1-x))


def finite_M(M,nu):
    assert M>=3 and nu>=0
    pre=[(c,n-c) for n in range(1,M) for c in range(n+1)]
    post=[(c,n-d-c,d) for n in range(1,M) for d in range(1,n+1) for c in range(n-d+1)]
    pidx={z:i for i,z in enumerate(pre)}
    qidx={z:i for i,z in enumerate(post)}
    Q=lil_matrix((len(pre),len(pre)))
    K=lil_matrix((len(post),len(post)))
    B=lil_matrix((len(pre),len(post)))
    t=np.zeros(len(post))
    for i,(c,l,d) in enumerate(post):
        size=c+l+d
        t[i]=THETA*d
        K[i,i]-=t[i]
        for k,z,b,mort in ((0,c,BC,DC),(1,l,BL,DL),(2,d,BD,DD)):
            if z==0: continue
            br=z*b; dr=z*mort
            K[i,i]-=(br+dr)
            if size<M-1:  # birth below establishment boundary
                nxt=[c,l,d];nxt[k]+=1
                K[i,qidx[tuple(nxt)]]+=br
            if size>1 and not (k==2 and d==1):
                nxt=[c,l,d];nxt[k]-=1
                K[i,qidx[tuple(nxt)]]+=dr
            # otherwise death of last D destroys commitment potential
        K[i,i]-=nu*l
        if l:
            K[i,qidx[(c+1,l-1,d)]]+=nu*l
    for i,(c,l) in enumerate(pre):
        size=c+l
        for k,z,b,mort in ((0,c,BC,DC),(1,l,BL,DL)):
            if z==0:continue
            br=z*b;dr=z*mort
            Q[i,i]-=br+dr
            if size<M-1:
                nxt=[c,l];nxt[k]+=1
                Q[i,pidx[tuple(nxt)]]+=br
                # mark newborn daughter as D; original parent remains
                B[i,qidx[(c,l,1)]]+=br
            if size>1:
                nxt=[c,l];nxt[k]-=1
                Q[i,pidx[tuple(nxt)]]+=dr
        Q[i,i]-=nu*l
        if l: Q[i,pidx[(c+1,l-1)]]+=nu*l
    q=spsolve((-K).tocsc(),t)
    A=spsolve((-Q).tocsc(),B.tocsr()@q)
    AC,AL=A[pidx[(1,0)]],A[pidx[(0,1)]]
    qdirect=q[qidx[(0,0,1)]]
    phi=founder_fraction(nu)
    D=2*((1-phi)*AC+phi*AL)/qdirect
    return D,AC,AL,qdirect


def finite_fast_repair_endpoint(M):
    _,AC,_,qd=finite_M(M,0.)
    return 2*AC/qd  # initial L repairs instantaneously


def unbounded_ratio(nu):
    FC=BC/(DC-BC)
    FL=(BL+nu*FC)/(DL-BL+nu)
    return 2*((1-founder_fraction(nu))*FC+founder_fraction(nu)*FL)


def run():
    FC=BC/(DC-BC)
    q_inf=1-(BD+DD+THETA-np.sqrt((BD+DD+THETA)**2-4*BD*DD))/(2*BD)
    assert np.isclose(FC,.6)
    assert np.isclose(q_inf,np.sqrt(4.25)-1.5)
    root_inf=brentq(lambda v:unbounded_ratio(v)-1,0,200,xtol=1.e-12)
    assert abs(root_inf-3.42421226182291)<1.e-9
    print('q_infinity:',format(q_inf,'.12f'))
    print('unbounded D(0), D(infty), crossing:',format(unbounded_ratio(0),'.12f'),1.2,format(root_inf,'.12f'))
    print('M | D(0) | D(fast repair) | crossover (if demonstrated) | qM')
    for M in (3,4,5,6,8,10):
        D0,AC,AL,qd=finite_M(M,0.)
        Dfast=finite_fast_repair_endpoint(M)
        if D0<1<Dfast:
            root=brentq(lambda v:finite_M(M,v)[0]-1,0,1000,xtol=1e-10)
        else:
            root=None
        print(f'{M} | {D0:.9f} | {Dfast:.9f} | {root} | {qd:.9f}')
        assert D0<Dfast
        assert 0<qd<=q_inf
        assert D0<=unbounded_ratio(0)
        assert Dfast<=2*FC+1.e-12
        assert finite_M(M,100000)[0] < Dfast+1e-5
        if root is not None: assert abs(finite_M(M,root)[0]-1)<1e-9
    # One example test of explicit exponential upper bound for post-hit direct success.
    eta=2.0
    assert eta*BC<DC and eta*BL<DL
    assert BD*(eta-1)+DD*(1/eta-1)<THETA
    for M in (3,4,5,6,8,10):
        qd=finite_M(M,0.)[3]
        assert q_inf-qd<=eta**(1-M)+1e-12
        a=min(DC-BC,DL-BL); Bmax=max(BC,BL); m=M//2
        delta=(Bmax/a)*(eta**(m-M)+m*(eta-1)/(eta**m-1))
        for nu in (0.,0.1,0.7,3.,25.,100.):
            _,AC,AL,_=finite_M(M,nu)
            FC=BC/(DC-BC)
            FL=(BL+nu*FC)/(DL-BL+nu)
            assert -1e-9<=q_inf*FC-AC<=delta+1e-9
            assert -1e-9<=q_inf*FL-AL<=delta+1e-9
    print('All generator, root, endpoint, direct-bound and insertion-bound checks passed.')

if __name__=='__main__':run()
