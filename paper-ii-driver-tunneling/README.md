# Paper II — Repair-dependent driver-pathway competition

Companion computational package for B. Hardisty, *Repair-Dependent Competition Between Direct Driver Founding and Stochastic Tunneling Before Clonal Establishment* (JMB manuscript; November 2026).

**Integration status:** This draft pull request currently contains the entry points, dependency list, one representative exact-expression check, and an independent CI workflow. The remaining checked source files, figure generators, TikZ schematic, and figure assets must be added from the local `Paper2_Companion_GitHub_Ready.zip` package **before merge**. The absence of these files is expected to cause CI to fail while this draft is incomplete.

When source import is complete, run from this directory:

```bash
python -m pip install -r requirements.txt
python run_checks.py --full
python make_figures.py
```

Scripts are separated from the Paper I project; no Paper I source file or workflow has been edited. Figure 1 is TikZ; the other six manuscript/appendix figures are reproduced from Python. Numerical curves illustrate behavior; rigorous finite-M witness claims rely on exact-rational residual checks, which certify inequalities at designated parameter values but not complete root counts. The manuscript takes the rare-driver coefficient at fixed M before considering M→∞.

**Publication note:** A Paper-II-specific release tag and archived DOI have **not** been created. Do not cite a DOI until a release/Zenodo deposit is verified.
