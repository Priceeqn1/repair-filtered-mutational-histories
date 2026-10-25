"""Exact rational check of killed-chain -> composition weights normalization."""
from fractions import Fraction as F
from math import comb,factorial
from itertools import combinations

r,g,e,nu=F(4,5),F(13,10),F(1,2),F(2,3)

for k in range(1,11):
    lam=lambda i:e*(g**i-r**i)
    pi=lambda i,j:F(comb(i,j))*((r/g)**j)*((1-r/g)**(i-j))/(1-(r/g)**i)
    s=lambda i:lam(i)/(lam(i)+i*nu)
    E=[F(1)]
    for i in range(1,k+1):
        E.append(s(i)*sum(pi(i,j)*E[j] for j in range(i)))
    a=lambda j:e*r**j/(lam(j)+j*nu)
    Z=F(0)
    for subset_size in range(k):
        for S in combinations(range(1,k),subset_size):
            bs=(0,)+S+(k,)
            weight=F(factorial(k))
            for t in range(len(bs)-1): weight/=factorial(bs[t+1]-bs[t])
            for j in S: weight*=a(j)
            Z+=weight
            path=F(1)
            # downward from k to previous levels
            for curr,prev in reversed(tuple(zip(bs[1:],bs[:-1]))):
                path*=s(curr)*pi(curr,prev)
            baseline=e*(g-r)**k/(lam(k)+k*nu)
            assert path==baseline*weight,('route identity',k,S)
    assert E[k]==e*(g-r)**k/(lam(k)+k*nu)*Z,('normalization',k)
print('PASS: exact Feynman-Kac/route normalization for every path up to k=10 (1023 routes), rational data')
