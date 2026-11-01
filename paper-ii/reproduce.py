#!/usr/bin/env python3
"""Paper II: verified founder-source calculations and auxiliary diagnostic figure."""
import argparse
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def pristine(nu):
    if not (math.isfinite(nu) and nu >= 0):
        raise ValueError("nu must be finite and nonnegative")
    # Stable positive solution of (2+4nu)x^2-x-4nu=0.
    return (1 + math.sqrt(1 + 16*nu*(2+4*nu)))/(2*(2+4*nu))

def inverse_nu(x):
    if not 0.5 <= x < 1:
        raise ValueError("physical x must be in [1/2,1)")
    return x*(2*x-1)/(4*(1-x*x))

def inverse_derivative(x):
    return (4*x-1-x*x)/(4*(1-x*x)**2)

def check():
    assert pristine(0) == 0.5
    for nu in [0, 1e-6, 0.01, 0.1, 1, 10, 100]:
        x = pristine(nu)
        assert 0.5 <= x < 1.0, (nu,x)
        assert math.isclose((2+4*nu)*x*x-x-4*nu,0,abs_tol=2e-12)
        assert math.isclose(inverse_nu(x),nu,rel_tol=2e-11,abs_tol=2e-12)
        assert 0.5 < 1/(1+x) <= 2/3
        h = min(1e-6, (1-x)/4)
        if h > 0:
            central = (inverse_nu(x+h)-inverse_nu(x-h))/(2*h)
            assert math.isclose(central,inverse_derivative(x),rel_tol=1e-5)
        assert inverse_derivative(x)>0
    print("PASS: Paper II founder-source algebra checks")

def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    (ROOT/"figures").mkdir(exist_ok=True)
    (ROOT/"data").mkdir(exist_ok=True)
    nus = [10**(-4+6*i/300)-1e-4 for i in range(301)]
    nus = [max(0.0,v) for v in nus]
    phis = [1/(1+pristine(v)) for v in nus]
    with (ROOT/"data"/"founder_source_diagnostic.csv").open("w",newline="") as fh:
        wr=csv.writer(fh)
        wr.writerow(["nu","x_pristine","phi_lesion_founder"])
        wr.writerows((v,pristine(v),p) for v,p in zip(nus,phis))
    fig,ax=plt.subplots(figsize=(6.5,4))
    ax.plot(nus,phis)
    ax.set(xlabel="Repair intensity nu",ylabel="Lesion-founder fraction phi",
           title="Auxiliary diagnostic: exact homeostatic founder source")
    ax.set_xscale("symlog",linthresh=0.01)
    fig.tight_layout()
    for ext in ("pdf","png"):
        fig.savefig(ROOT/"figures"/("founder_source_diagnostic."+ext),dpi=180)
    plt.close(fig)
    print("PASS: auxiliary diagnostic figure and CSV generated")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group(required=True)
    for option in ("checks","figures","all"):
        group.add_argument("--"+option,action="store_true")
    args=ap.parse_args()
    if args.checks or args.all: check()
    if args.figures or args.all: figures()
