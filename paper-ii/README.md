# Paper II: repair-dependent stochastic tunneling

Independent reproducibility workspace for *Repair-Dependent Competition Between Direct Driver Founding and Stochastic Tunneling Before Clonal Establishment* (manuscript v11).

This folder is deliberately separate from the Paper I code and workflow. **Current coverage is preliminary:** it verifies the exact one-locus founder-source quadratic, its inverse coordinate and derivative, and illustrates the founder composition. It does **not** yet reconstruct the full finite-state killed generators, sequential/direct ratios, manuscript figures, or the M=9/M=10 numerical claims. Those require porting and checking the original figure calculations before treating the repository as full manuscript reproduction.

## Run

Python 3.10+:

```bash
python paper-ii/reproduce.py --checks
python -m pip install -r paper-ii/requirements.txt
python paper-ii/reproduce.py --all
```

`--checks` uses only Python's standard library. `--figures` produces an **auxiliary diagnostic** PDF and PNG in `paper-ii/figures/`, and an accompanying CSV in `paper-ii/data/`. It is not represented as an article figure.

## Manuscript formulas being checked

For the illustrative one-locus parameters `(U,r,g,epsilon)=(1/2,1/2,2,1/4)`, Eq. (67) in manuscript v11 reads

```
(2+4 nu) x^2 - x - 4 nu = 0,   1/2 <= x < 1
phi = 1/(1+x)
```

The inverse coordinate in Eq. (68) is `nu=x(2x-1)/(4(1-x*x))`; its derivative is `(4x-1-x*x)/(4(1-x*x)^2)`. The checks use an independent quadratic root and inverse calculation.

## Pending full reproduction

- Build and independently cross-check the transient pre-hit/post-hit finite-state generators at each M.
- Reproduce exact `A_C,A_L,q_C,q_L`, the founder-weighted ratio `D_M=2S/Z`, and the manuscript's finite-M plots.
- Verify the M=10 values reported in manuscript §9.1 (approximately 1.007919, 0.962283, 1.016642 at nu = 0, 0.1, 10) **from the generator**, not by merely comparing hard-coded numbers.
- Compare M=9 with M=10 without asserting unverified M=9 stationary points.
- Record the exact figure-to-script and parameter mapping, then add regression fixtures for every displayed figure.

See `docs/REPRODUCIBILITY_SCOPE.md` for coverage and limitations.
