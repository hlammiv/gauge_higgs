#!/usr/bin/env python3
"""Harvest the L_s=20 wedge run and resolve the q>=5 Coulomb wedge with a ROBUST photon-mass estimator.

The driver's m2 is a linear R^-1 = a + b*phat2 fit over the 4 lowest shells, which was degenerate at L_s=16.
Here we REFIT from the emitted R(p) shells with two cross-checks:
  - massless DIAGNOSTIC (model-free): in Coulomb R(p) ~ Z/phat2  => R*phat2 = const (flat in p).
    In Higgs R ~ Z/(phat2+m2) => R*phat2 RISES with p. We report the ratio
        S = (R*phat2)[2nd shell] / (R*phat2)[1st shell].  S~1 massless (Coulomb); S>>1 massive (Higgs).
  - m2 from the lowest TWO shells: b=(y2-y1)/(x2-x1), a=y1-b*x1, m2=a/b   (y=1/R, x=phat2).
A point is MASSLESS (Coulomb) if the 2-pt m2 is consistent with 0 (|m2|<m2_lo) AND S<S_thr; MASSIVE (Higgs)
if m2_2pt>m2_hi and S>S_thr; else AMBIGUOUS. Wedge = massless cells inside the matter-condensed (<cos> high)
region for q>=5.  q=2 control should go MASSIVE at large kappa (no wedge)."""
import glob, re, os
import numpy as np

H   = re.compile(r"Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
M2  = re.compile(r"m2=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
COS = re.compile(r"<cos>_link=([\-\d.eE]+)")
SG1 = re.compile(r"sigma1=([\-\d.eEna]+)")
SHELL = re.compile(r"phat2=([\d.]+)\s+R=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

def load(dirs):
    recs = {}
    for d in dirs:
        for fn in glob.glob(d + "/*.out"):
            t = open(fn).read(); h = H.search(t)
            if not h: continue
            shells = [(f(a), f(b), f(c)) for a, b, c in SHELL.findall(t)]
            shells = [s for s in shells if np.isfinite(s[1]) and s[1] > 0]
            shells.sort()
            m = M2.search(t); c = COS.search(t); s1 = SG1.search(t)
            recs[(int(h[5]), float(h[3]), float(h[4]))] = dict(    # (q, beta, kappa)
                Ls=int(h[1]), Lt=int(h[2]),
                m2_drv=f(m[1]) if m else float('nan'),
                cos=f(c[1]) if c else float('nan'),
                s1=f(s1[1]) if s1 else float('nan'),
                shells=shells)
    return recs

def robust(r):
    sh = r['shells']
    if len(sh) < 2: return dict(m2=float('nan'), S=float('nan'))
    (x1, R1, _), (x2, R2, _) = sh[0], sh[1]
    y1, y2 = 1.0/R1, 1.0/R2
    b = (y2 - y1) / (x2 - x1)
    a = y1 - b * x1
    m2 = a / b if abs(b) > 1e-12 else float('nan')
    S = (R2 * x2) / (R1 * x1) if R1 * x1 != 0 else float('nan')   # massless diagnostic ~1
    return dict(m2=m2, S=S)

M2_LO, M2_HI, S_THR, COS_THR = 0.02, 0.05, 1.5, 0.5
def state(r):
    rb = robust(r); m2, S = rb['m2'], rb['S']
    if not np.isfinite(m2) or not np.isfinite(S): return "amb", rb
    massless = (abs(m2) < M2_LO) and (S < S_THR)
    massive  = (m2 > M2_HI) and (S > S_THR)
    return ("massless" if massless else "massive" if massive else "amb"), rb

if __name__ == "__main__":
    recs = load(("u1f_wedge20", "u1f_wedge20_lenore"))
    print(f"# L_s=20 wedge points: {len(recs)}")
    qs = sorted(set(k[0] for k in recs)); betas = sorted(set(k[1] for k in recs)); kaps = sorted(set(k[2] for k in recs))
    print(f"# q={qs} beta={betas} kappa={kaps}")
    print("\n## photon state vs kappa (deconfined side). massless=Coulomb wedge ; massive=Higgs ; q=2=no-wedge control")
    for q in qs:
        print(f"\n q={q}:")
        for b in betas:
            cells = []
            for k in kaps:
                r = recs.get((q, b, k))
                if not r: cells.append(f"k{k:g}:----"); continue
                st, rb = state(r)
                tag = {"massless":"M0", "massive":"Mx", "amb":"??"}[st]
                cells.append(f"k{k:g}:{tag}(m2={rb['m2']:+.2f},S={rb['S']:.1f})")
            print(f"   b={b}: " + "  ".join(cells))
    # wedge verdict: fraction massless at kappa>=1.5 in the matter-condensed region, per q
    print("\n## WEDGE VERDICT: massless fraction at kappa>=1.5 & <cos> condensed (Higgs region):")
    for q in qs:
        hi = [recs[(q,b,k)] for b in betas for k in kaps
              if (q,b,k) in recs and k >= 1.5 and recs[(q,b,k)]['cos'] > COS_THR]
        nm = sum(1 for r in hi if state(r)[0] == "massless")
        print(f"   q={q}: {nm}/{len(hi)} massless deep in the Higgs region"
              + ("   <-- Coulomb WEDGE" if hi and nm/len(hi) > 0.4 else ""))
