#!/usr/bin/env python3
"""Standalone exact counterexample certificates for Paper I.
All claims use Fraction arithmetic; decimal/float comparison is avoided.
Run: python verify_referee_certs.py
"""
from fractions import Fraction as F
from math import comb, factorial


def weights(k, a):
    c1 = sum(comb(k, j) * a[j] for j in range(1, k))
    c2 = sum(F(factorial(k), factorial(i)*factorial(j-i)*factorial(k-j))
             * a[i] * a[j] for i in range(1, k) for j in range(i+1, k))
    return c1, c2

def gap_singleton(k, center, a):
    c1,c2=weights(k,a)
    lo=F(comb(k,center))*a[center]/c1
    up=sum(F(factorial(k),factorial(min(center,j))*factorial(abs(center-j))*factorial(k-max(center,j)))
           * a[center]*a[j] for j in range(1,k) if j!=center)/c2
    return lo-up

def partition(k, r, g, eps, nu):
    a=[F(1)]+[eps*r**j/(eps*(g**j-r**j)+j*nu) for j in range(1,k)]
    Z=[F(0)]*(k+1)
    for i in range(1,k+1):
        Z[i]=F(1)+sum(comb(i,j)*a[j]*Z[j] for j in range(1,i))
    return a, Z[k]

def main():
    uniform = {j:F(1) for j in range(1,12)}
    assert gap_singleton(12,6,uniform)==F(1750,350037)

    for nu in (F(0),F(1)):
        a={j:F(2**j,4**j-2**j+2*j*nu*3**j) for j in range(1,14)}
        gap=gap_singleton(14,4,a)
        assert gap>0
        if nu==0:
            assert gap==F(1017173459848427695564485247810291575,
                          180754611939941975863563191643324144974)
        print('Hall exact gap sign for nu=',nu,':',gap>0,
              '(numerator digits:',len(str(gap.numerator)),')')

    k=15; r=F(9,10); g=F(5); eps=F(1,2)
    tails={}
    for nu, low, high in [(F(0),134850,134851),(F(1,5),134758,134759)]:
        a,Z=partition(k,r,g,eps,nu)
        t=(F(1)+F(k,2)*(a[1]+a[k-1]))/Z
        # All inequalities below are exact integer comparisons.
        assert F(low,10**6)<t<F(high,10**6)
        tails[nu]=t
        print('event tail rational enclosure for nu=',nu,': ',low,'/10^6 < T <',high,'/10^6')
    assert tails[F(1,5)]<tails[F(0)]
    print('ALL CERTIFICATES PASSED (no floating-point arithmetic).')

if __name__=='__main__':
    main()
