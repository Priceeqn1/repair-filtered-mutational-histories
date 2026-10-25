"""Independent exact-integer/rational validation of the new revision results.
Run: python checks/verify_new_theorems.py. Python standard library only.
This corroborates but is not a proof of the infinite-family theorems.
"""
from fractions import Fraction as F
from itertools import combinations
from math import factorial, comb
from collections import defaultdict


def mul(A,B):
    o=defaultdict(int)
    for i,a in A.items():
        for j,b in B.items(): o[tuple(x+y for x,y in zip(i,j))]+=a*b
    return dict(o)

def sub(A,B):
    c=defaultdict(int,A)
    for k,v in B.items():c[k]-=v
    return {k:v for k,v in c.items() if v}

checks=0
for k in range(2,8):
    dim=k-1
    C={};T={}
    for m in range(1,k+1):
        c=defaultdict(int)
        tails={(q,l):defaultdict(int) for q in range(1,k) for l in range(1,k+1)}
        for S in combinations(range(1,k),m-1):
            b=(0,)+S+(k,)
            cf=factorial(k)
            for u,v in zip(b[:-1],b[1:]):cf//=factorial(v-u)
            mon=tuple(int(j in S) for j in range(1,k))
            c[mon]+=cf
            for q in range(1,k):
                L=sum(j<=q for j in S)
                for l in range(1,L+1):tails[q,l][mon]+=cf
        C[m]=dict(c)
        for a,v in tails.items(): T[m,*a]=dict(v)
    for m in range(1,k):
        for q in range(1,k):
            for l in range(1,m+1):
                delta=sub(mul(C[m],T[m+1,q,l]),mul(C[m+1],T[m,q,l]))
                assert all(v>=0 for v in delta.values()),('prefix_tail',k,m,q,l)
                checks+=1
print('PASS coefficientwise prefix-tail inequalities:',checks,'through k=7')

n=5;r=F(4,5);g=F(13,10);eps=F(1,2)
for U in [F(1,100),F(1,8),F(3,5)]:
    for nu in [F(0),F(2,5),F(3)]:
        def M(i,j):return F(comb(n-i,j-i))*U**(j-i)*(1-U)**(n-j)
        def D(j):return eps*g**j+j*nu-eps*r**j*(1-U)**(n-j)
        def a(j):return eps*r**j*(1-U)**(n-j)/D(j)
        for k in range(1,n+1):
            for m in range(1,k+1):
                for S in combinations(range(1,k),m-1):
                    path=(0,)+S+(k,)
                    exact=F(1)
                    for u,v in zip(path[:-1],path[1:]):exact*=eps*r**u*M(u,v)/D(v)
                    gaps=[v-u for u,v in zip(path[:-1],path[1:])]
                    c=F(factorial(k))
                    for d in gaps:c/=factorial(d)
                    rhs=U**k*comb(n,k)*eps*(1-U)**(n-k)/D(k)*c
                    for j in S:rhs*=a(j)
                    assert exact==rhs,('finite_U',U,nu,k,path)
print('PASS exact finite-U route factorization, 3 U values x 3 repair values, all routes through k=5')
