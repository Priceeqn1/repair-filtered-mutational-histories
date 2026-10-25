#!/usr/bin/env python3
"""Exact rational check of the new finite-k score counterexample."""
from fractions import Fraction as F
from itertools import combinations
from math import factorial

k=4; r=F(1,7); g=F(12,7); eps=F(1,2); nu=F(0)
h={j:F(j,1)/(eps*(g**j-r**j)+j*nu) for j in range(1,k)}
a={j:eps*r**j/(eps*(g**j-r**j)+j*nu) for j in range(1,k)}
def c(S):
    cuts=(0,)+S+(k,)
    v=factorial(k)
    for p,q in zip(cuts,cuts[1:]): v//=1  # preserve an integer before rational division
    denom=1
    for p,q in zip(cuts,cuts[1:]): denom*=factorial(q-p)
    return F(v,denom)
def ave(m):
    mass=F(0); weighted=F(0)
    for S in combinations(range(1,k),m-1):
        w=c(S)
        for j in S: w*=a[j]
        mass+=w; weighted+=w*sum((h[j] for j in S),F(0))
    return weighted/mass
Delta=ave(4)-ave(3)
assert Delta == F(252,209), Delta
assert h[1] == F(14,11),h[1]
assert Delta/h[1] == F(18,19)
assert r+g<2
print('PASS: exact k=4 counterexample; Delta/h1 = 18/19 < 1; Delta>0.')
