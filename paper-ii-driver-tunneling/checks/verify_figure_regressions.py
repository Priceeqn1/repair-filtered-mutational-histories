"""Numerical-data regressions for Paper II manuscript figures (v11).

Fixed benchmarks transcribed from the manuscript and independent certification logs.
Compare recalculated mathematical values, not PDF bytes or styling.
"""
import math
import struct
import sys
from pathlib import Path
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"reproduce"))
sys.path.insert(0,str(ROOT/"reproduce"/"sensitivity"))

from verify_sign_classification import finite_ratio, unbounded
from reproduce_sensitivity import finite_ratio as sensitivity_ratio
from verify_generator import model, pop

def near(actual, expected, tol, label):
    assert math.isfinite(float(actual)) and abs(float(actual)-expected)<=tol, (
        f"{label}: actual={actual}, expected={expected}, tolerance={tol}")
    print(f"PASS {label}: {float(actual):.10g}")

def figure_2():
    near(brentq(lambda v:pop(v)-1,0,1),.192345334861229,1e-9,"Fig 2 population crossing")
    near(brentq(lambda v:model(v)[0]-1,0,3.5),1.10964042855,1e-8,"Fig 2 conditional crossing")

def figure_3():
    # Independently compare the five-state model with manuscript's exact
    # nu=0 opportunity and establishment fractions (M=3; birth rates 10,4).
    import numpy as np
    from scipy.linalg import solve
    states=[(0,1),(1,0),(0,2),(1,1),(2,0)]
    idx={z:i for i,z in enumerate(states)}
    G=np.zeros((5,5));eligible=np.zeros(5);establish=np.zeros(5)
    for i,z in enumerate(states):
        size=sum(z)
        for k,b in enumerate((10.,4.)):
            if z[k]==0:continue
            rate=z[k]*b
            G[i,i]-=rate
            if size==2: establish[i]+=rate
            else:
                zz=list(z);zz[k]+=1;G[i,idx[tuple(zz)]]+=rate
                eligible[i]+=rate
            rate=z[k]
            G[i,i]-=rate
            if size>1:
                zz=list(z);zz[k]-=1;G[i,idx[tuple(zz)]]+=rate
    near(solve(-G,eligible)[idx[(0,1)]],20/21,1e-11,"Fig 3 eligible births at nu=0")
    near(solve(-G,establish)[idx[(0,1)]],16/21,1e-11,"Fig 3 establishment probability at nu=0")

def figures_4_5():
    # Source: verified M9/M10 finite-generator rational enclosures.
    for M,nu,val in [
      (9,0,.998471800785),(9,.1,.964530797011),(9,6,1.009204963908),
      (10,0,1.007918564),(10,.1,.962282795),(10,10,1.016641534),
      (10,100,1.009430346)]:
        near(finite_ratio(M,nu),val,8e-6,f"Figs 4-5 M={M} nu={nu}")
    # The M=9 curve crosses parity again by nu=1000, verified separately
    # by exact rational residuals in certify_M10_parity.py.
    assert finite_ratio(9,1000)<1

def figure_6():
    for U,theta,nu,expected in [
      (.2,.001,0,1.006962849),
      (.2,.001,.1,.995593758),
      (.2,.001,10,1.011302765),
      (.1,.01,0,1.004034951),
      (.1,.01,.2,.999464301),
      (.1,.01,10,1.009540599)]:
        near(sensitivity_ratio(10,nu,U=U,theta_L=theta),
             expected,2e-6,f"Fig 6 U={U} theta={theta} nu={nu}")

def output_files():
    names=[
      "population_repair_switch_M3",
      "opportunity_tradeoff_M3",
      "reentrant_pathways",
      "threshold_regimes",
      "finite_M10_sensitivity",
      "lesion_input_only",
    ]
    for name in names:
        pdf=ROOT/"figures"/(name+".pdf")
        png=ROOT/"figures"/(name+".png")
        assert pdf.exists() and pdf.stat().st_size>1500, str(pdf)
        assert pdf.open("rb").read(5)==b"%PDF-",str(pdf)
        assert png.exists() and png.stat().st_size>1500,str(png)
        with png.open("rb") as f:
            header=f.read(24)
        assert header[:8]==b"\x89PNG\r\n\x1a\n"
        width,height=struct.unpack(">II",header[16:24])
        assert width>=300 and height>=200,(name,width,height)
        print(f"PASS generated PDF/PNG: {name} ({width} x {height})")
    # Figure 1 is manuscript-integrated TikZ, not made by Python.
    tex=ROOT/"figures"/"model_schematic.tex"
    assert tex.exists() and tex.stat().st_size>200

if __name__=="__main__":
    figure_2();figure_3();figures_4_5();figure_6();output_files()
    print("PASS: all v11 figure numerical regressions and file checks")
