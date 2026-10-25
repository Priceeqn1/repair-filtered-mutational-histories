#!/usr/bin/env python3
"""Exact arithmetic certificates for weighted refinement Hall inequalities.

No dependencies beyond the Python standard library.
For k<=6 each polynomial Hall slack has nonnegative integer coefficients
for arbitrary positive activities. Therefore an inclusion-preserving coupling
exists for every pair of adjacent conditional route-rank measures.

Also verifies rational counterexamples for k=12 (uniform activities) and for
k=14 with activities arising from the homeostatic repair model.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb, factorial


def strata(k, rank):
    """Pairs (bit-mask, exact multinomial coefficient), rank=|S|."""
    result = []
    for ss in combinations(range(1, k), rank):
        sizes = [b-a for a,b in zip((0,)+ss, ss+(k,))]
        c = factorial(k)
        for d in sizes:
            c //= factorial(d) # valid: divide successive factorials
        mask = sum(1<<(i-1) for i in ss)
        result.append((mask,c))
    return result


def hall_polynomial(k, rank, selected_indices):
    """Coefficient dictionary for C_low * w(upper-neighbors(A)) - C_high * w(A).

    Monomials represented by (union-mask, intersection-mask): bit in
    intersection means exponent 2; union-only bits mean exponent 1.
    """
    low, high = strata(k,rank), strata(k,rank+1)
    A = set(selected_indices)
    neigh = {j for j,(t,_) in enumerate(high)
             if any((low[i][0]&t)==low[i][0] for i in A)}
    polynomial=defaultdict(int)
    for s,c in low:
        for j in neigh:
            t,d=high[j]
            polynomial[(s|t,s&t)] += c*d
    for i in A:
        s,c=low[i]
        for t,d in high:
            polynomial[(s|t,s&t)] -= c*d
    return {mon:coef for mon,coef in polynomial.items() if coef}


def verify_through(max_k=6):
    total=0
    for k in range(2,max_k+1):
        for rank in range(1,k-1):
            n = len(strata(k,rank))
            examined = 0
            # skip empty/full families (trivial equalities)
            for bitset in range(1,(1<<n)-1):
                A = [i for i in range(n) if (bitset>>i)&1]
                poly = hall_polynomial(k,rank,A)
                negative = [(mon,v) for mon,v in poly.items() if v<0]
                if negative:
                    raise AssertionError((k,rank,A,negative[:3]))
                examined+=1
            total+=examined
            print(f"k={k}, |S|={rank}->{rank+1}: {examined} exact Hall polynomials OK")
    print(f"PASS: {total} nontrivial exact polynomial inequalities certified for k<= {max_k}.")


def route_pair_coefficient(k,i,j):
    if i>j: i,j = j,i
    return factorial(k)//(factorial(i)*factorial(j-i)*factorial(k-j))


def singleton_hall_difference(k,center,a):
    """P_{|S|=1}(S={center}) - P_{|S|=2}(center in S)."""
    C1=sum(comb(k,i)*a[i] for i in range(1,k))
    C2=sum(route_pair_coefficient(k,i,j)*a[i]*a[j]
           for i in range(1,k) for j in range(i+1,k))
    lower=comb(k,center)*a[center]/C1
    upper=sum(route_pair_coefficient(k,center,i)*a[center]*a[i]
              for i in range(1,k) if i!=center)/C2
    return lower,upper,lower-upper


def counterexamples():
    k=12
    a={i:Fraction(1) for i in range(1,k)}
    lo,hi,gap=singleton_hall_difference(k,6,a)
    assert lo==Fraction(462,2047)
    assert hi==Fraction(868,3933)
    assert gap==Fraction(1750,350037) and gap>0
    print(f"k=12, all a_j=1: P(rank1 contains 6)={lo}, P(rank2 contains 6)={hi}, gap={gap}.")

    k=14
    for nu in (0,1):
        # r=2/3, g=4/3, eps=1/2 (repair rate nu).
        a={i:Fraction(2**i,4**i-2**i+2*i*nu*3**i)
           for i in range(1,k)}
        lo,hi,gap=singleton_hall_difference(k,4,a)
        assert gap>0
        print(f"k=14, r=2/3,g=4/3,eps=1/2,nu={nu}: "
              f"P_2(V_4)= {float(lo):.12f}, P_3(V_4)= {float(hi):.12f}, "
              f"exact rational gap > 0: {gap>0} ({float(gap):.12f}).")

    # Failure of *coefficientwise proof method*, not failure of Hall inequality.
    low = strata(7,1)
    center_index=[i for i,(mask,c) in enumerate(low) if mask==(1<<(3-1))][0]
    polynomial=hall_polynomial(7,1,[center_index])
    key = ((1<<0)|(1<<2)|(1<<3),0) # a1*a3*a4
    assert polynomial[key]==-245
    print("k=7, |S|=1->2, A={{3}}: coefficient of a1*a3*a4 is -245; "
          "coefficientwise strategy stops working (not a counterexample).")


if __name__=='__main__':
    verify_through(6)
    counterexamples()
