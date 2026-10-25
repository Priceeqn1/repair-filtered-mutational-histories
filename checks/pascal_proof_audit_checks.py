"""Independent exact audit of Pascal TP2 and prefix ordering.
Uses direct enumeration of break-point subsets, not the manuscript recurrence.
Pure Python; integer coefficients in polynomial ring Z[a_1,...,a_{k-1}].
"""
from itertools import combinations
from math import factorial
from collections import defaultdict


def poly_add(A, B):
    C = defaultdict(int, A)
    for mon, c in B.items(): C[mon] += c
    return {x:v for x,v in C.items() if v}


def poly_mul(A,B):
    C = defaultdict(int)
    for a,aa in A.items():
        for b,bb in B.items():
            C[tuple(x+y for x,y in zip(a,b))] += aa*bb
    return {x:v for x,v in C.items() if v}


def poly_sub(A,B):
    return poly_add(A,{m:-v for m,v in B.items()})


def route_polys(k):
    """Q_{k,m} and U_{k,m,q}, direct independent route enumeration."""
    n = k-1
    C, U, M = {}, {}, {}
    for m in range(1,k+1):
        Cp = defaultdict(int)
        Up = [defaultdict(int) for _ in range(k-1)]
        Mp = [[defaultdict(int) for _ in range(k)] for _ in range(m)]
        for S in combinations(range(1,k),m-1):
            bs = (0,)+S+(k,)
            gaps = [bs[j+1]-bs[j] for j in range(m)]
            c = factorial(k)
            for g in gaps: c //= factorial(g)  # WRONG if divided sequentially? exact each step yes
            exp = tuple(int(i in S) for i in range(1,k))
            Cp[exp] += c
            for q in range(1,k):
                Up[q-1][exp] += c * sum(i<=q for i in S)
            for r,j in enumerate(S,1):
                Mp[r][j][exp] += c
        C[m] = dict(Cp)
        for q in range(1,k): U[m,q]= {mon:v for mon,v in Up[q-1].items() if v}
        for r in range(1,m):
            for j in range(1,k):M[m,r,j]={mon:v for mon,v in Mp[r][j].items() if v}
    return C,U,M

# First verify the multinomial integer formula used above independently.
for k in range(2,10):
    for m in range(1,k+1):
        for S in combinations(range(1,k),m-1):
            bs=(0,)+S+(k,)
            den=1
            for j in range(m):den*=factorial(bs[j+1]-bs[j])
            a=factorial(k)//den
            b=factorial(k)
            for j in range(m): b//=factorial(bs[j+1]-bs[j])
            assert a==b

count_det = count_prefix = count_break = 0
for k in range(2,9):
    C,U,M=route_polys(k)
    zero=tuple(0 for _ in range(k-1))
    # Note: Q_{n,s} for target n uses same activity ring with padded coordinates.
    qs={n:route_polys(n)[0] for n in range(1,k+1)}
    def pad(P,n): return {mon+ (0,)*(k-n):v for mon,v in P.items()}
    for u in range(1,k):
        for v in range(u+1,k+1):
            for s in range(1,k+1):
                Pu = pad(qs[u].get(s,{}),u)
                Pus=pad(qs[u].get(s+1,{}),u)
                Pv = pad(qs[v].get(s,{}),v)
                Pvs=pad(qs[v].get(s+1,{}),v)
                det=poly_sub(poly_mul(Pu,Pvs),poly_mul(Pv,Pus))
                assert all(val>=0 for val in det.values()),('TP2',k,u,v,s,det)
                count_det+=1
    for m in range(1,k):
        for q in range(1,k):
            D=poly_sub(poly_mul(C[m],U[m+1,q]),poly_mul(C[m+1],U[m,q]))
            assert all(val>=0 for val in D.values()),('prefix',k,m,q,D)
            count_prefix+=1
        for r in range(1,m):
            for q in range(1,k):
                Va= {}
                Vb={}
                for j in range(1,q+1):
                    Va=poly_add(Va, M[m,r,j]);Vb=poly_add(Vb,M[m+1,r,j])
                D=poly_sub(poly_mul(C[m],Vb),poly_mul(C[m+1],Va))
                assert all(val>=0 for val in D.values()),('breakpoint',k,m,r,q,D)
                count_break+=1
    print('exact direct-enumeration checks passed through k=',k, flush=True)
print('PASS: TP2 determinants',count_det,'prefix polynomial inequalities',count_prefix,'breakpoint lower-tail polynomial inequalities',count_break)
# Additional independent check of the ordered-breakpoint factorization.
count_factor=0
for k in range(2,10):
    C,U,M=route_polys(k)
    for m in range(2,k+1):
        for r in range(1,m):
            for j in range(1,k):
                if (j<r or k-j<m-r):
                    assert not M[m,r,j]
                    continue
                prefix=route_polys(j)[0][r]
                # reverse local suffix activity b_t=a_{k-t} for t<k-j.
                suffix=route_polys(k-j)[0][m-r]
                P={}
                for pexp,pcoeff in prefix.items():
                    for sexp,scoeff in suffix.items():
                        e=[0]*(k-1)
                        for i,v in enumerate(pexp,1):e[i-1]+=v
                        for i,v in enumerate(sexp,1):e[k-i-1]+=v
                        e[j-1]+=1
                        e=tuple(e)
                        P[e]=P.get(e,0)+ (factorial(k)//factorial(j)//factorial(k-j))*pcoeff*scoeff
                assert P==M[m,r,j],('factorization',k,m,r,j,P,M[m,r,j])
                count_factor+=1
print('PASS: breakpoint factorization for', count_factor, 'nonempty cases through k=9')
