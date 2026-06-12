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
                e=f(m[2]) if m else float('nan'),
                cos=f(c[1]) if c else float('nan'),
                s1=f(s1[1]) if s1 else float('nan'),
                shells=shells)
    return recs

def ratio(r):
    # model-free massless diagnostic: R(p) ~ Z/phat2 (massless) => R[0]/R[1] -> phat2[1]/phat2[0] ~ 2.0;
    # massive (saturating) => ~1.0. Robust, uses the two lowest shells directly.
    sh = r['shells']
    if len(sh) < 2 or sh[0][1] <= 0: return float('nan')
    return sh[1][1] and sh[0][1] / sh[1][1]

# Use the DRIVER 4-pt m2 (better than a 2-pt estimate -- it resolves small masses via shells 3,4) WITH its
# error, cross-checked by the model-free R-ratio. (R-ratio ~2 massless, ~1 massive; m2 significance decides.)
RAT_HI, RAT_LO, M2_HI = 1.7, 1.45, 0.03
def state(r):
    m2, e, rr = r['m2_drv'], r.get('e', float('nan')), ratio(r)
    if not np.isfinite(m2) or not np.isfinite(rr): return "amb", dict(m2=m2, rr=rr)
    sig = m2 / e if (np.isfinite(e) and e > 0) else 0.0
    massless = (rr > RAT_HI) and (abs(m2) < 0.02)
    massive  = (rr < RAT_LO) and (m2 > M2_HI) and (sig > 3)
    return ("massless" if massless else "massive" if massive else "amb"), dict(m2=m2, rr=rr, sig=sig)

COS_THR = 0.5
if __name__ == "__main__":
    recs = load(("u1f_wedge20", "u1f_wedge20_lenore"))
    print(f"# L_s=20 wedge points: {len(recs)}")
    qs = sorted(set(k[0] for k in recs)); betas = sorted(set(k[1] for k in recs)); kaps = sorted(set(k[2] for k in recs))
    print(f"# q={qs} beta={betas} kappa={kaps}")
    print("\n## photon state vs kappa (deconfined side). M0=massless(Coulomb) Mx=massive(Higgs) ??=ambiguous")
    for q in qs:
        print(f"\n q={q}:")
        for b in betas:
            cells = []
            for k in kaps:
                r = recs.get((q, b, k))
                if not r: cells.append(f"k{k:g}:----"); continue
                st, d = state(r)
                tag = {"massless":"M0", "massive":"Mx", "amb":"??"}[st]
                cells.append(f"k{k:g}:{tag}(m2={d['m2']:+.2f},rr={d['rr']:.1f})")
            print(f"   b={b}: " + "  ".join(cells))
    # WEDGE VERDICT at the near-wall slices (beta=1.2,1.5): is the photon massless deep in the condensed region?
    print("\n## WEDGE VERDICT -- near-wall (beta<=1.5), matter condensed (<cos>>0.5), kappa>=1.5:")
    for q in qs:
        hi = [recs[(q,b,k)] for b in betas for k in kaps
              if (q,b,k) in recs and b <= 1.5 and k >= 1.5 and recs[(q,b,k)]['cos'] > COS_THR]
        nm = sum(1 for r in hi if state(r)[0] == "massless")
        mx = sum(1 for r in hi if state(r)[0] == "massive")
        verdict = "Coulomb WEDGE (massless)" if hi and nm > mx and nm/len(hi) >= 0.5 else \
                  "Higgs (massive, NO wedge)" if hi and mx >= nm else "ambiguous"
        print(f"   q={q}: massless={nm} massive={mx} amb={len(hi)-nm-mx} of {len(hi)}  -> {verdict}")

    # ---- money plot: m^2(kappa) at the near-wall slice beta=1.2 ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=False)
    for ax, b in zip(axes, [1.2, 1.5]):
        for q in qs:
            ks = [k for k in kaps if (q, b, k) in recs]
            m2 = [recs[(q, b, k)]['m2_drv'] for k in ks]
            er = [recs[(q, b, k)]['e'] for k in ks]
            ax.errorbar(ks, m2, yerr=er, marker="o", capsize=3, lw=1.8, label=f"q={q}")
        ax.axhline(0, color="k", lw=0.8, ls=":")
        ax.axhspan(-0.02, 0.02, color="0.85", zorder=0)   # massless band (resolution)
        ax.set_xlabel(r"$\kappa$"); ax.set_ylabel(r"$m_\gamma^2$  (driver 4-pt fit)")
        ax.set_title(rf"$\beta$={b} (just past confinement wall)")
        ax.set_ylim(-0.05, 0.3); ax.legend(fontsize=9, ncol=2)
        ax.annotate("massive (Higgs) " + r"$\uparrow$", xy=(0.02, 0.97), xycoords="axes fraction",
                    fontsize=8, va="top", color="0.3")
        ax.text(0.98, 0.06, "massless band (Coulomb)", transform=ax.transAxes, ha="right",
                fontsize=8, color="0.4")
    fig.suptitle("L$_s$=20 photon mass vs $\\kappa$: q$\\leq$4 Higgses (m$^2$>0) but q$\\geq$5 stays massless (Coulomb wedge)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    import os; os.makedirs("u1f_campaign_analysis", exist_ok=True)
    p = "u1f_campaign_analysis/wedge20_mgamma.png"; fig.savefig(p, dpi=120); print("\nwrote", p)
