"""Exact-rational audit: frozen-source branching first moments vs rare-U route weights.

Run: python verify_lineage_correspondence.py
Only Python standard library is required.
"""
from fractions import Fraction as F
from math import comb, factorial
from itertools import combinations

n=5; r=F(4,5); g=F(13,10); eps=F(1,2); nu=F(2,5)

def mutation(i,j,U):
    if j<i: return F(0)
    return F(comb(n-i,j-i))*U**(j-i)*(1-U)**(n-j)

def route_paths(k):
    return [(0,)+tuple(S)+(k,) for m in range(k) for S in combinations(range(1,k),m)]

def no_repair_route_count(path,U):
    res=F(1)
    for u,v in zip(path[:-1],path[1:]):
        D=eps*g**v+v*nu-eps*r**v*mutation(v,v,U)
        res*=eps*r**u*mutation(u,v,U)/D
    return res

def direct_weights(k):
    weights={}
    for path in route_paths(k):
        S=path[1:-1]
        lens=[b-a for a,b in zip(path[:-1],path[1:])]
        c=F(factorial(k))
        for d in lens: c/=factorial(d)
        for j in S: c*=eps*r**j/(eps*(g**j-r**j)+j*nu)
        weights[path]=c
    return weights

def full_means(U):
    # Solve A z = -b by exact rational Gaussian elimination.
    a=[]; b=[]
    for i in range(1,n+1):
        row=[]
        for j in range(1,n+1):
            v=F(0)
            if j<=i: v+=eps*r**j*mutation(j,i,U)
            if j==i: v-=eps*g**i+i*nu
            if j==i+1: v+=(i+1)*nu
            row.append(v)
        a.append(row)
        b.append(-eps*mutation(0,i,U))
    for x in range(n):
        pivot=a[x][x]
        assert pivot
        for j in range(x+1,n):
            ratio=a[j][x]/pivot
            for l in range(x,n):a[j][l]-=ratio*a[x][l]
            b[j]-=ratio*b[x]
    z=[F(0)]*n
    for x in reversed(range(n)):
        z[x]=(b[x]-sum((a[x][l]*z[l] for l in range(x+1,n)),F(0)))/a[x][x]
    return [F(1)]+z

def run():
    k=3
    weights=direct_weights(k)
    target=sum(weights.values(),F(0))
    lam=eps*(g**k-r**k)
    leading=eps*comb(n,k)/(lam+k*nu)*target
    print('k=',k,'route count=',len(weights),'leading abundance coefficient=',float(leading))
    prev_cond_err = None
    prev_repair_share = None
    prev_route_errors = {}
    for U in [F(1,50),F(1,100),F(1,200),F(1,400)]:
        all_k=full_means(U)[k]
        paths={path:no_repair_route_count(path,U) for path in weights}
        nr=sum(paths.values(),F(0))
        assert nr<=all_k, 'repair-free type count must not exceed total'
        assert all_k>0
        direct=paths[(0,k)]/all_k
        direct_limit=weights[(0,k)]/target
        cond_err=abs(direct-direct_limit)
        repair_share=(all_k-nr)/all_k
        print(f'U={float(U):.5f}, E[Zk]/N/U^k={float(all_k/U**k):.9f}, '
              f'repair-containing share={float(repair_share):.7f}, '
              f'conditional direct={float(direct):.9f}, '
              f'limit={float(direct_limit):.9f}, error={float(cond_err):.3g}')
        if prev_cond_err is not None:
            assert cond_err < prev_cond_err, 'route law must converge as U decreases'
            assert repair_share < prev_repair_share, 'repair-containing share must decline'
        prev_cond_err = cond_err
        prev_repair_share = repair_share
        for path, w in weights.items():
            route_leading=eps*comb(n,k)/(lam+k*nu)*w
            route_err=abs(paths[path]/U**k-route_leading)
            if path in prev_route_errors:
                assert route_err < prev_route_errors[path], f'route {path} did not converge'
            prev_route_errors[path]=route_err
    assert abs(full_means(F(1,400))[k]/F(1,400)**k-leading)<F(3,20)
    print('PASS exact-rational mean-field route checks')

if __name__=='__main__':run()
