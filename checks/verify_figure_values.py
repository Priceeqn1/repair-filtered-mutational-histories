"""Regression tests for the three manuscript figure calculations.

Only the Python standard library is required. Tests numerical claims, not
rendering/pixel output and not the general mathematical theorems.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from repair_histories.model import (two_hit_stats, event_distribution,
                                     survival_factor, equal_leading_burden_mutation_ratio)
from math import isclose

def approx(x,y,tol=1e-8):
    assert isclose(x,y,rel_tol=tol,abs_tol=tol),(x,y)

def main():
    p0,b0=two_hit_stats(0)
    p1,b1=two_hit_stats(1)
    approx(p0,5/21)
    approx(p1,25/41)
    approx(100*b1/b0,8.118811881188119)
    for k in (2,3,7,10):
        approx(survival_factor(k,0),1,1e-10)
        for nu in (0,1):
            p=event_distribution(k,nu)
            approx(sum(p),1,1e-12)
            assert all(x>0 for x in p)
    p0=event_distribution(10,0)
    p1=event_distribution(10,1)
    mean0=sum((i+1)*p for i,p in enumerate(p0))
    mean1=sum((i+1)*p for i,p in enumerate(p1))
    approx(mean0,4.5648,5e-5)
    approx(mean1,3.2732,5e-5)
    # Likelihood ratio p1(m)/p0(m) decreases with event count.
    ratios=[a/b for a,b in zip(p1,p0)]
    assert all(ratios[i]>ratios[i+1] for i in range(9))
    e7=survival_factor(7,1)
    approx(e7,0.00956703934,1e-8)
    rr=equal_leading_burden_mutation_ratio(7,0,1)
    approx(rr,1.94294429,1e-7)
    print('PASS figure regressions:',
          f'two-hit burden={100*b1/b0:.8f}%',
          f'two-hit direct shares={100*two_hit_stats(0)[0]:.5f}%/{100*two_hit_stats(1)[0]:.5f}%',
          f'k10 means={mean0:.6f},{mean1:.6f}',
          f'E7(1)={e7:.11f}',f'U2/U1={rr:.8f}')

if __name__=='__main__':main()
