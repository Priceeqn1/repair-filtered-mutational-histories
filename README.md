# Repair-filtered mutational histories: reproducibility code

**Companion code for** B. Hardisty, *Repair Filters Mutational Histories by Event Count: A Pascal Ordering Theorem with a Lineage Interpretation* (research manuscript, 2026; prepared for the *Bulletin of Mathematical Biology*).

This is a **standalone project**, deliberately separate from the author's earlier Haldane-related manuscript repositories. It reproduces the three main manuscript figures and runs independent finite exact-arithmetic checks of the combinatorial, Feynman–Kac, and genealogical calculations.

> **Scope.** The stochastic route distributions describe leading rare-mutation equilibrium contributions. They acquire a reproductive-ancestry interpretation only in the specified frozen-background branching–immigration model. The code does not reconstruct observed cancer pedigrees, estimate cancer incidence, or replace the manuscript's general mathematical proofs.

## Reproduce everything

Python **3.10+** recommended. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate               # Windows: .venv\\Scripts\\activate
python -m pip install -r requirements-figures.txt
python reproduce.py --all
```

The exact mathematical checks run using only the standard library, **without installing NumPy or Matplotlib**:

```bash
python reproduce.py --checks
```

The figures can be regenerated separately:

```bash
python reproduce.py --figures
```

Scripts can be invoked from **any working directory**, because paths are resolved relative to their own files. Successful completion prints `PASS`. Individual checks and figure scripts are also directly runnable with `python path/to/script.py`.

## Figures generated

| Manuscript figure | File in `figures/` | Source script | Calculation |
|---|---|---|---|
| Fig. 1, three-panel overview | `paper_I_three_panel_accurate.pdf` | `scripts/generate_three_panel_figure.py` | Direct/sequential diagram; two-hit bars; exact normalized 10-hit event-count probabilities |
| Fig. 2, two-hit mechanism | `paper_I_two_hit_mechanism.pdf` | `scripts/generate_result_figures.py` | Leading total two-hit abundance and conditional direct-route share versus repair |
| Fig. 3, equal-burden contrasting histories | `paper_I_jump_count.pdf` | `scripts/generate_result_figures.py` | Normalized event-count probabilities for `k=7` at `nu=0` and `nu=1` |

All figures are exported as **vector PDFs** and PNG previews. Their plotted values are written to `data/figure_*.csv` so readers can check the numbers without extracting graphics. Styling need not be pixel-identical across Matplotlib versions; the equations and data are the reproducibility targets.

Manuscript parameters for all three figures: `r = 0.8`, `g = 1.3`, `eps = 0.5`; repair levels `nu = 0, 1` where compared, `k = 10` in Fig. 1C and `k = 7` in Fig. 3.

The files use the **same figure basenames** as the LaTeX manuscript, making it possible to replace the corresponding PDFs directly in the Overleaf `figures/` directory if desired.

## Exact checks

The scripts in `checks/` test:

- Coefficientwise Pascal/TP2, breakpoint, prefix and suffix inequalities for finite instances;
- Exact rational killed-chain / route normalization;
- Exact marked-lineage, finite-`U` and rare-mutation calculations;
- Hall-type route-coupling certificates and a failure of set inclusion;
- The event-size stochastic-order counterexample;
- Numeric regression values quoted in the three figures.

See [`docs/REPRODUCIBILITY_SCOPE.md`](docs/REPRODUCIBILITY_SCOPE.md) for scope and [`docs/FIGURE_MAPPING.md`](docs/FIGURE_MAPPING.md) for equations and expected numbers.

## Repository layout

```text
repair-filtered-mutational-histories/
├── README.md
├── reproduce.py
├── requirements-figures.txt
├── src/repair_histories/          # shared finite route/repair recursions
├── scripts/                       # generators for all three figures
├── figures/                       # generated vector PDFs and PNGs
├── data/                          # plotted CSV data and local check log
├── checks/                        # exact-arithmetic verification scripts
├── docs/                          # model, figure and publication notes
└── .github/workflows/             # reproducibility CI
```

## Citation and provenance

The current repository accompanies **Paper I**, not the related preprint *Repair-Mediated Path Filtering of Mutation Load in Multilocus Homeostatic Populations* (Hardisty, 2026, DOI [`10.64898/2026.10.11.758278`](https://doi.org/10.64898/2026.10.11.758278)). That preprint develops the original equilibrium/repair-filtering recursion; Paper I studies weighted route ordering and genealogical interpretations.

See [`CITATION.cff`](CITATION.cff). When a public repository and archived release exist, add the actual repository URL and a release DOI, if minted. **Do not cite the parent preprint DOI as the code repository DOI.**

## Before publishing

Choose and add a `LICENSE` file appropriate to how you want others to reuse this code. No license is asserted here without the author's approval. The GitHub URL and any Zenodo DOI must be filled in only after a repository/release actually exists.

For manual publishing, see [`docs/GITHUB_PUBLISHING.md`](docs/GITHUB_PUBLISHING.md).
