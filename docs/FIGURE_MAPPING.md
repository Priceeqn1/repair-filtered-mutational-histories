# Figure formulas and expected values

All three figures derive from

\[
\lambda_j=\varepsilon(g^j-r^j),\quad
 a_j(\nu)=\frac{\varepsilon r^j}{\lambda_j+j\nu},\quad
 C_{n,1}=1,\quad C_{n,m}=\sum_{j=1}^{n-1}\binom nj a_j(\nu) C_{j,m-1},
\]

with invalid or out-of-range `C` interpreted as zero. The conditional count law is `P(J=m|k,nu)=C[k,m]/sum_s C[k,s]`. We have **not** sampled a stochastic ancestry simulation to make these figures; they are exact finite-recursion calculations in floating-point arithmetic.

## Figure 1: overview

Panel A is a conceptual cartoon (0→2 and 0→1→2). Panel B displays the two-hit leading-order ratios described below. Panel C shows the normalized event-count distribution for `k=10` under `nu=0` and `nu=1` using `r=.8, g=1.3, eps=.5`. Conditional means are approximately `4.565` and `3.273`, respectively. The plot does not claim to show exact finite-`U` reconstructed tumor lineages.

## Figure 2: two-hit abundance and direct route share

For `lambda_i = eps*(g^i-r^i)`, suppressing the common factor `binom(n,2) U^2` in the leading two-hit abundance gives

\[
  B_2^{[2]}(\nu)=\frac{\varepsilon}{\lambda_2+2\nu}
   \left(1+\frac{2\varepsilon r}{\lambda_1+\nu}\right),\qquad
  P_{\mathrm{direct}}(\nu)=\frac{\lambda_1+\nu}{\lambda_1+\nu+2\varepsilon r}.
\]

With the stated parameters, `B_2^{[2]}(1)/B_2^{[2]}(0)=0.08118811881188119`, `P_direct(0)=5/21=0.238095238...`, and `P_direct(1)=25/41=0.609756097...`.

## Figure 3: equal leading abundance, different route laws

`k=7`, same `r,g,eps`, `nu_1=0`, `nu_2=1`. The repair-filter factor is

\[
 E_k(\nu)=\frac{\varepsilon(g-r)^k}{\lambda_k+k\nu}
\sum_{m=1}^k C_{k,m}(\nu).
\]

Values: `E_7(0)=1`, `E_7(1)≈0.00956703934`, and the mutation-rate ratio matching **leading** `k`-hit abundance is

\[
 \frac{U_2}{U_1}=\left[\frac{E_7(0)}{E_7(1)}\right]^{1/7}
 \approx 1.94294429.
\]

The manuscript separately proves that nearby finite mutation rates can be chosen to match the **exact** stationary burden via an implicit-function argument. Figure 3 itself compares the limiting conditional route laws; it does **not** solve the full finite-`U` equilibrium matching problem numerically.
