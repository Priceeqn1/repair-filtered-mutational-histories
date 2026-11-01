# Auditable rational-residual certificate for finite-`M` pathway parity

**Companion to:** `reproduce/certify_M10_parity.py`. This note explains precisely what its rational arithmetic certifies. It concerns the specific two-lesion-state, two-driver-type branching process encoded in that script, at the listed rational repair rates and thresholds `M=9,10`. It does **not** count all zeros of the pathway ratio or prove that every parameter choice exhibits reversals.

## 1. Transient generators and exact targets

The transient pre-hit state space `P_M` consists of nonzero vectors `(n_C,n_L)` with total `n_C+n_L<M`. The post-hit state space `H_M` consists of nonzero vectors `(n_C,n_L,m_C,m_L)` with total size less than `M` and `m_C+m_L>0`. A jump that leaves a transient state space is an absorbing exit and retains its rate in the diagonal loss term, even if it does not appear as an off-diagonal transition. The generator construction uses the **finite, killed** process: births, deaths, repair, and post-hit irreversible commitment.

Write `Q` and `K` for the transient generators on `P_M` and `H_M`, respectively, `B` for the nonnegative pre-to-post second-hit insertion-rate matrix, and `t` for the post-hit commitment-rate vector. Define

```text
(-K) q = t,                 (-Q) A = B q .
```

For a pristine and a lesion-bearing singleton founder, the entries `q_C,q_L` of `q`, and `A_C,A_L` of `A`, are selected by the corresponding singleton indices. In the homeostatic one-locus source, `x` is the positive root of `(2+4ν)x² - x - 4ν = 0`. The certified sequential/direct ratio is

```text
D_M(ν) = 2 (x A_C + A_L) / (x q_C + q_L).
```

All rates defining `Q,K,B,t` and the tested repair rates `ν` are rational; numerical linear solves provide only *candidate* solutions.

## 2. A post-hit bound: a uniform commitment hazard

Every transient post-hit state contains at least one two-driver cell. Its irreversible commitment intensity is therefore at least `θ_min = 1/1,000,000` (the smaller of the two two-driver commitment intensities). With `1` the all-ones vector,

```text
(-K) 1 >= θ_min 1.
```

Since `-K` is an invertible transient killed-generator M-matrix, `(-K)^{-1}` is entrywise nonnegative. Multiplying this inequality by that inverse yields

```text
0 <= (-K)^{-1} 1 <= (1/θ_min) 1.
```

Thus `||(-K)^{-1}||∞ <= 1/θ_min`. If `q̃` is any rational vector and its **exactly evaluated** rational residual is `r_q=t+K q̃`, then

```text
||q - q̃||∞ <= ||r_q||∞ / θ_min =: E_q.
```

Note that no bound on the original floating-point solver error is assumed: the solver is only used to obtain `q̃`, whose rational residual is calculated afresh.

## 3. A pre-hit bound: the population-size barrier

Both one-driver cell types are strictly subcritical. Let

```text
a = min(d_C-b_C, d_L-b_L) > 0.
```

For any pre-hit state `z`, set `V(z)=|z|/a`, with `|z|=n_C+n_L`, so `1 <= |z| <= M-1`. Repair preserves population size. Away from the upper establishment boundary, the generator of total population size is `Σ_k (b_k-d_k)n_k <= -a |z|`. A birth that exits the transient set at total size `M` contributes `-b_k n_k V(z)` to `Q V` rather than a larger positive increment, so the killing boundary can only *strengthen* the inequality. Deaths to extinction likewise use `V(0)=0`. Therefore,

```text
(-Q) V >= 1   componentwise on P_M.
```

By positivity of the killed Green operator,

```text
0 <= (-Q)^{-1} 1 <= V <= ((M-1)/a) 1,
||(-Q)^{-1}||∞ <= (M-1)/a.
```

Let `Ã` be a rational approximate solution generated with the candidate `q̃`. Define its exact rational residual `r_A = B q̃ + Q Ã`, and let

```text
||B||∞ = max_z Σ_y B(z,y) .
```

The error obeys

```text
(-Q)(A-Ã) = B(q-q̃) + r_A
||A-Ã||∞ <= ((M-1)/a) ( ||r_A||∞ + ||B||∞ E_q ) =: E_A .
```

This remains valid when the computed candidates are inaccurate: the exact residual will then make the enclosure correspondingly wide.

## 4. Enclose the physical homeostatic root

At `ν=0`, `x=1/2` exactly. For every other tested nonnegative rational `ν`, the code computes a floating-point *guess* and constructs rational endpoints `x_-<x_+` near it. It then evaluates the polynomial

```text
p_ν(x)=(2+4ν)x²-x-4ν
```

in exact rational arithmetic, checking

```text
1/2 < x_- < x_+ < 1,  p_ν(x_-) < 0 < p_ν(x_+).
```

On this interval `p_ν` is strictly increasing, so the unique physical root satisfies `x_-<x<x_+`. The float is *not* used to certify this inclusion; the endpoint signs are.

## 5. Enclose the pathway ratio and certify strict signs

Select `ã_C,ã_L,q̃_C,q̃_L` from the rational candidate vectors. The componentwise errors from Sections 2–3 imply, for both founder types,

```text
ã_i-E_A <= A_i <= ã_i+E_A,
q̃_i-E_q <= q_i <= q̃_i+E_q.
```

For the endpoint formula below, also verify `ã_C-E_A>0` and `q̃_C-E_q>0`; these are explicit guard assertions in the script. The source interval then gives

```text
N_- = 2 [ x_- (ã_C-E_A) + (ã_L-E_A) ],
N_+ = 2 [ x_+ (ã_C+E_A) + (ã_L+E_A) ],
Z_- = x_- (q̃_C-E_q) + (q̃_L-E_q),
Z_+ = x_+ (q̃_C+E_q) + (q̃_L+E_q).
```

When `N_->0` and `Z_->0`, the resulting *exact rational* certificate is

```text
N_-/Z_+ <= D_M(ν) <= N_+/Z_- .
```

The code decides `D_M>1` only when the **lower** endpoint is above one, and `D_M<1` only when the **upper** endpoint is below one. Displayed twelve-digit decimals are for readability; the sign tests themselves compare rational numbers.

Strict signs at finitely many repair values plus continuity of the finite-state absorption probabilities establish **at least two crossings** in the cited examples, but not their locations, uniqueness, or total number.

## 6. Exact-scope limitations and reproducibility

* The state-space sizes and generator convention in `prepare(M,ν)` are part of the certificate's hypotheses. A change in the threshold boundary, second-hit insertion convention, or commitment rule requires a new argument and rerun.
* Floating-point arithmetic chooses candidate linear-system solutions and root approximations. It does not decide the certified parity signs; rational residuals, root-bracket sign tests, and rational interval endpoints do.
* Invertibility/positivity uses killed finite-state Markov generators. The strictly subcritical one-driver types justify the pre-hit bound; positive post-hit commitment justifies the post-hit bound.
* The result is a finite list of **strict sign certificates**, not an exhaustive polynomial root-isolation result, general parameter theorem, or independent formal proof that the coded model exactly matches the manuscript.

Run from the repository root:

```bash
python paper-ii-driver-tunneling/reproduce/certify_M10_parity.py
python paper-ii-driver-tunneling/checks/verify_all_eight_published.py
```

The second script checks the eight stated manuscript enclosures against the separately reported bounds.
