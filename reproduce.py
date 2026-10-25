#!/usr/bin/env python3
"""Paper I reproducibility entrypoint (run from any current directory).

Usage:
    python reproduce.py --all       # figures plus exact checks
    python reproduce.py --figures   # requires matplotlib and numpy
    python reproduce.py --checks    # Python standard library only
"""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
FIGURES=['scripts/generate_three_panel_figure.py','scripts/generate_result_figures.py']
CHECKS=sorted(str(p.relative_to(ROOT)) for p in (ROOT/'checks').glob('*.py'))

def run(files):
    for f in files:
        print('\n===',f,'===',flush=True)
        subprocess.run([sys.executable,str(ROOT/f)],cwd=ROOT,check=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--all',action='store_true')
    g.add_argument('--figures',action='store_true')
    g.add_argument('--checks',action='store_true')
    args=ap.parse_args()
    if args.figures or args.all: run(FIGURES)
    if args.checks or args.all: run(CHECKS)
    print('\nPASS: selected reproduction steps completed.')

if __name__=='__main__': main()
