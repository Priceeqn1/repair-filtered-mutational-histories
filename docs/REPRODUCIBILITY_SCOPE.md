# Reproducibility scope and limitations

## What is computationally verified

The provided checks are finite-domain, automated tests. They use Python's arbitrary-precision integers and/or exact `fractions.Fraction` arithmetic for statements described as exact. The numerical figure generators use double-precision floating-point arithmetic, with normalization and regression assertions. The scripts document their checked ranges.

## What remains mathematical, rather than computational, proof

The manuscript's claims for **all** target sizes, arbitrary positive activities, or complete continuous parameter regions rely on the written induction and stochastic-ordering proofs. No finite test suite can establish them universally.

## Genealogical interpretation

The decreasing auxiliary chain is a Feynman–Kac representation, not a chronological cell lineage. The manuscript constructs a separate, subcritical frozen-background branching–immigration realization. Its rare-mutation marked ancestry coefficient agrees with the weighted-route coefficient. No claim is made about empirical tumor phylogenies, selection of irreversible oncogenic drivers, cancer incidence, or universal stochastic realizations of the homeostatic ODE.

## Figures

Figure 1A is schematic. Figures 1B, 1C, 2, and 3 use the explicit analytical recursions under displayed fixed parameters. PDFs are regenerated using Matplotlib; font rasterization and PDF metadata may change across operating systems or library versions. The generated CSVs are the preferred numerical reference.
