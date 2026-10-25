"""Dependency-free mathematical calculations shared by manuscript figures.

Definitions (all parameters dimensionless):
  lambda_i = eps*(g**i - r**i)
  a_i(nu) = eps*r**i / (lambda_i + i*nu)
  C[n,m] = delta(m,1) + sum_{j=1}^{n-1} binom(n,j)*a_j*C[j,m-1]
  P(J=m | k,nu) = C[k,m]/sum_m C[k,m].

The final quantity represents weighted route contributions conditional
on target k. It is not an empirical tumor-event distribution.
"""
from math import comb

DEFAULT_R = 0.8
DEFAULT_G = 1.3
DEFAULT_EPS = 0.5


def check_parameters(k, nu, r, g, eps):
    if not isinstance(k, int) or k < 1:
        raise ValueError('k must be a positive integer')
    if not (0 < r < 1 < g and eps > 0 and nu >= 0):
        raise ValueError('Require 0<r<1<g, eps>0, nu>=0')


def lambda_i(i, *, r=DEFAULT_R, g=DEFAULT_G, eps=DEFAULT_EPS):
    if i < 1: raise ValueError('i must be positive')
    return eps * (g**i - r**i)


def route_activity(i, nu, *, r=DEFAULT_R, g=DEFAULT_G, eps=DEFAULT_EPS):
    check_parameters(i, nu, r, g, eps)
    return eps * r**i/(lambda_i(i,r=r,g=g,eps=eps)+i*nu)


def route_partition(k, nu, *, r=DEFAULT_R, g=DEFAULT_G, eps=DEFAULT_EPS):
    """Return unnormalized C[k,m] for m=1,..,k.

    Dynamic program: O(k^3) time, O(k^2) memory, no path enumeration.
    """
    check_parameters(k, nu, r, g, eps)
    C = [[0.0]*(k+1) for _ in range(k+1)]
    for n in range(1, k+1):
        C[n][1] = 1.0
        for m in range(2, n+1):
            C[n][m] = sum(comb(n,j)*route_activity(j,nu,r=r,g=g,eps=eps)*C[j][m-1]
                          for j in range(m-1,n))
    return C[k][1:k+1]


def event_distribution(k, nu, *, r=DEFAULT_R, g=DEFAULT_G, eps=DEFAULT_EPS):
    """Normalized limiting event-count law P(J=m|k), m=1,...,k."""
    weights = route_partition(k,nu,r=r,g=g,eps=eps)
    total = sum(weights)
    return [x/total for x in weights]


def survival_factor(k, nu, *, r=DEFAULT_R, g=DEFAULT_G, eps=DEFAULT_EPS):
    """E_k(nu), the repair-filter factor relative to nu=0."""
    weights = route_partition(k,nu,r=r,g=g,eps=eps)
    return eps*(g-r)**k*sum(weights)/(lambda_i(k,r=r,g=g,eps=eps)+k*nu)


def two_hit_stats(nu, *, r=DEFAULT_R, g=DEFAULT_G, eps=DEFAULT_EPS):
    """Return (direct-route share, leading burden coefficient / binom(n,2)).

    The ratio of the second return value at nu=1 vs nu=0 is the manuscript's
    normalized total two-hit abundance, to leading order in mutation U.
    """
    check_parameters(2, nu, r, g, eps)
    lam1 = lambda_i(1,r=r,g=g,eps=eps)
    lam2 = lambda_i(2,r=r,g=g,eps=eps)
    direct_share = (lam1+nu)/(lam1+nu+2*eps*r)
    burden_coef = eps/(lam2+2*nu)*(1+2*eps*r/(lam1+nu))
    return direct_share, burden_coef


def equal_leading_burden_mutation_ratio(k, nu1, nu2, *, r=DEFAULT_R,g=DEFAULT_G,eps=DEFAULT_EPS):
    """Leading U2/U1 for equal k-hit burdens, fixed target k."""
    if nu2 < nu1: raise ValueError('Require nu2>=nu1')
    return (survival_factor(k,nu1,r=r,g=g,eps=eps)/
            survival_factor(k,nu2,r=r,g=g,eps=eps))**(1.0/k)
