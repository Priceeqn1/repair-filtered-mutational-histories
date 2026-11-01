# Paper II integration — staging note

This repository's existing root workflow and files remain **Paper I**.

Paper II's proposed companion directory is `paper-ii-driver-tunneling/`, with its own Python requirements, verification commands and GitHub Actions workflow.

## Import remaining source files before merging this draft PR

The locally verified package is `Paper2_Companion_GitHub_Ready.zip` (previously shared in ChatGPT). Its `paper-ii-driver-tunneling/` directory contains the complete Python source tree, TikZ schematic, figures, and documentation. Add the **contents** to this branch, preserving relative paths. Do not upload only the ZIP if the goal is browsable/reviewable source on GitHub.

This draft currently includes the verified launcher files and a representative exact-expression check, not the full source suite. The workflow is expected to fail until the remaining source files have been added.

## After source import

Run `python paper-ii-driver-tunneling/run_checks.py --full` and `python paper-ii-driver-tunneling/make_figures.py` and review the GitHub Actions run before marking this PR ready or merging. The smooth numerical sweeps are not complete algebraic root certifications.

## Publication

Create a Paper-II-specific release tag only after reproducibility review. Confirm which Zenodo version DOI is produced; the Paper-I archive/DOI must not be described as a Paper-II archive.
