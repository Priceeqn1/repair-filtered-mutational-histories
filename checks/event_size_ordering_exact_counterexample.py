"""Exact counterexample to full stochastic jump-size ordering.

For the rare-mutation route weights, choose a route with probability
proportional to its composition weight, then select one of its jumps
uniformly. All arithmetic below uses Python's exact Fraction class.

Admissible model parameters: k=15, r=9/10, g=5, epsilon=1/2.
They satisfy 0<r<1<g, g+r>2, 0<epsilon<1.

Run: python event_size_ordering_exact_counterexample.py
No third-party packages are required.
"""
from fractions import Fraction as Q
from math import comb

K = 15
R = Q(9, 10)
G = Q(5)
EPS = Q(1, 2)


def calculate(nu):
    """Exact partition Z, exposure B, and P(D_event >= k-1)."""
    a = [Q(1)] + [
        EPS * R**j / (EPS * (G**j - R**j) + j * nu)
        for j in range(1, K)
    ]
    h = [Q(0)] + [
        Q(j) / (EPS * (G**j - R**j) + j * nu)
        for j in range(1, K)
    ]

    # Last-jump decomposition. Z_i contains all route weights to i,
    # normalized by the direct route's weight (which is 1).
    Z = [Q(0)] * (K + 1)
    B = [Q(0)] * (K + 1)  # B_i = -dZ_i / dnu
    for i in range(1, K + 1):
        Z[i] = Q(1) + sum(comb(i, j) * a[j] * Z[j]
                          for j in range(1, i))
        B[i] = sum(comb(i, j) * a[j] *
                   (B[j] + h[j] * Z[j]) for j in range(1, i))

    # Only the direct route and the two splits (1,k-1), (k-1,1)
    # have an event of size at least k-1.
    numerator = Q(1) + Q(K, 2) * (a[1] + a[K-1])
    minus_numerator_derivative = Q(K, 2) * (
        a[1] * h[1] + a[K-1] * h[K-1]
    )
    tail = numerator / Z[K]
    derivative = (numerator * B[K]
                  - minus_numerator_derivative * Z[K]) / Z[K]**2
    direct = Q(1) / Z[K]
    return tail, derivative, direct, Z[K]


if __name__ == "__main__":
    before = calculate(Q(0))
    after = calculate(Q(1, 5))
    # Strict inequalities are established with integers/rationals,
    # not by floating-point approximations.
    assert before[1] < 0
    assert after[0] < before[0]
    assert after[2] > before[2]
    assert G + R > 2
    for nu, result in [(Q(0), before), (Q(1, 5), after)]:
        tail, deriv, direct, z = result
        print(f"nu={nu}: P(D_event>=14)={float(tail):.12f}, "
              f"derivative={float(deriv):+.12f}, "
              f"P(D_event=15)={float(direct):.12f}, "
              f"Z={float(z):.12f}")
    print("Exact P(0.2)-P(0) is negative:", after[0] - before[0] < 0)
    print("All exact-rational assertions passed.")
