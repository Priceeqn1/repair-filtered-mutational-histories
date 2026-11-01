# Paper II — Repair-dependent driver-pathway competition

Companion computational and figure-generation code for:

**Repair-Dependent Competition Between Direct Driver Founding and Stochastic Tunneling Before Clonal Establishment** — Ben Hardisty (JMB manuscript; pre-submission working package, November 2026).

This folder is intended to live inside the **existing Paper I repository** [`repair-filtered-mutational-histories`](https://github.com/Priceeqn1/repair-filtered-mutational-histories), alongside—not in place of—its original scripts. Paper I studies repair-filtered mutational histories; Paper II examines competing direct and sequential driver-progression routes before a finite clone-size boundary. Neither paper's reproduction pipeline depends on the other's.

## Setup / run

From the `paper-ii-driver-tunneling` folder, using Python 3.10+:

```bash
python -m pip install -r requirements.txt
python run_checks.py              # fast numerical and algebraic checks
python run_checks.py --full       # add exact symbolic/rational certificates
python make_figures.py           # regenerate PDF+PNG figures in figures/
```

For paper-II-only scripts, no external data or network downloads are needed after Python dependencies are installed. `make_figures.py` runs the original plotting implementations plus a generator-based recreation of the exact population-switching plot.

## Figure-to-code map (JMB manuscript v11)

| Figure in manuscript | Generated file / source | Script | Method |
|---|---|---|---|
| Fig. 1: model and competing exits | `figures/model_schematic.tex` | LaTeX `\input{model_schematic.tex}` | **TikZ** (not Python) |
| Fig. 2: homeostatic population-weighted and conditional switching | `figures/population_repair_switch_M3.pdf` | `reproduce/plot_switches.py` | M=3 killed-generator inversion; exact root verified |
| Fig. 3: opportunity versus establishment | `figures/opportunity_tradeoff_M3.pdf` | `reproduce/plot_opportunity.py` | Five-state finite generator |
| Fig. 4: re-entrant pathway curves | `figures/reentrant_pathways.pdf` | `reproduce/plot_reentrant.py` | Unbounded analytic ratio + numerical finite generators |
| Fig. 5: finite-size/threshold regimes | `figures/threshold_regimes.pdf` | `reproduce/plot_threshold_regimes.py` | Finite-state generator sweeps |
| Fig. 6: biological sensitivity, M=10 | `figures/finite_M10_sensitivity.pdf` | `reproduce/plot_biological_sensitivity.py` | Finite-generator sweeps; marked rates verified separately |
| Appendix sensitivity figure | `figures/lesion_input_only.pdf` | `reproduce/plot_biological_sensitivity.py` | Unbounded analytic expression |
| Archived supplementary conditional figure (not in current PDF) | `figures/repair_switch_M3.pdf` | `reproduce/plot_switches.py` | M=3 generator |

`make_figures.py` generates numerically equivalent curves; plotting design may differ from the original manuscript PDF. For exact production artwork, use the figure PDFs in the matching Overleaf manuscript version. The TikZ schematic is sourced in its own `.tex` file.

## Verification mapping

- `reproduce/verify_generator.py` — direct 3-type generator, conditional formula, population crossover.
- `reproduce/verify_population_switch.py` — population-weighted exact expression.
- `reproduce/verify_balance.py` — exact opportunity–turnover–establishment identity and example.
- `reproduce/verify_common_green.py` — three-channel resolvent/sensitivity decomposition checks.
- `reproduce/verify_large_M_independent.py` — independent subcritical comparisons and finite-boundary calculations.
- `reproduce/verify_sign_classification.py` — unbounded derivative/sign expressions and finite-M numeric comparisons.
- `reproduce/verify_heterogeneous.py` and `reproduce/verify_finite_M.py` — independent post-hit lesion heterogeneity check; the original verification script had a temporary-file import which has been replaced with the displayed source-equilibrium formula.
- `reproduce/certify_M3_global.py` — exact polynomial/Descartes classification for specified M=3 four-type example.
- `checks/certify_M3_with_bracket.py` — additional exact-rational bracket check.
- `reproduce/certify_M10_parity.py` — rational-residual bounds for selected M=9 and M=10 rates.
- `checks/verify_all_eight_published.py` — independently checks all eight published M=9/10 witness intervals; tests the *printed* bounds, not just signs.
- `reproduce/sensitivity/certify_moderate.py`, `certify_lower.py` — strict rational-residual witnesses for sensitivity examples.

## Evidence and limits

1. **Analytic verification:** Some scripts verify closed-form identities symbolically or evaluate them independently against killed-generator inversions. They do not independently re-prove every theorem in the manuscript.
2. **Numerical curves:** Smooth finite-M sweeps are numerical evaluations. Their shapes, crossings, and inferred extrema are not globally certified by the plotted grid.
3. **Exact certificates:** Selected parity witnesses are enclosed by exact-rational residual inequalities after floating candidate solves. These certify strict inequalities at the designated points. In each finite-M example, sign changes and continuity force *at least* the asserted crossings; they do **not** count all real roots.
4. **M=3 classification:** For a specified rational four-type parameter example, the symbolic Descartes calculation establishes **no parity crossing and exactly one local maximum**, not universal classification for all M or rates.
5. **Large M:** The approximation theorem relies on its stated uniform subcriticality, positive commitment, and no-new-lesion assumptions; plots alone do not certify the uniform bounds.
6. **Biology:** All rate parameters are model illustrations, not clinical/tissue calibration. Establishment is an operational boundary, and commitment is a model-defined event.
7. **Mutation order:** Extract rare driver-mutation coefficients with fixed finite M before sending M to infinity. These scripts do not establish uniformity in the driver mutation parameter.

The figure source and verification scripts were gathered from the actual Paper II working-project archives and adjusted to run from this self-contained subdirectory; their underlying parameter families were not changed.

## Citations and permanent archiving

Publish this folder in the same repository as Paper I, but give **Paper II its own release tag and Zenodo *version DOI*** (if archiving through the linked Zenodo repository). If Zenodo connects a repository-wide concept DOI, carefully distinguish that concept DOI from a Paper-II-specific version record; don't claim the DOI is paper-specific until confirmed. Add the exact archived revision identifier to the paper's Code availability statement once released.

The repository owner should review code, generated outputs and license before making this public. The manuscript and repository should point to one another only when their URLs are verified.
