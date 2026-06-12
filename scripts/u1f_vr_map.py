#!/usr/bin/env python3
"""V(R)-based wedge map from the `pot` sweep -- the INDEPENDENT (position-space) companion to the
structure-factor wedge. For every (q,beta,kappa) it reconstructs V(R), Cornell-fits c+sigma*R-alpha/R, and
labels the phase by SHAPE: Confined sigma>0 (linear) | Coulomb sigma~0,alpha>0 (-1/R, charge-1 UNSCREENED) |
Higgs sigma~0,alpha~0 (FLAT, charge-1 SCREENED). Prints a per-q table + the wedge verdict, and draws V(R)
small-multiples (one panel per q). Reuses scripts/u1f_vr.py.
usage: scripts/u1f_vr_map.py u1f_vr/job_*.out"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import u1f_vr as vr

SIG_THR, ALPHA_THR = vr.SIG_THR, vr.ALPHA_THR
def phase(f):                       # 0 Confined, 1 Coulomb/wedge, 2 Higgs, -1 fail
    if f is None: return -1
    if f['sigma'] > SIG_THR:   return 0
    if f['alpha'] > ALPHA_THR: return 1
    return 2
NAME = {0: "Confined", 1: "Coulomb/wedge", 2: "Higgs(screened)", -1: "fit-fail"}

def main(files):
    recs = {}   # (q,beta,kappa) -> dict
    for fn in files:
        meta, W = vr.load(fn)
        Rs, VR, VRe = vr.V_of_R(W, meta['Rmax'])
        f = vr.fit(Rs, VR, VRe)
        recs[(meta['q'], meta['beta'], meta['kappa'])] = dict(meta=meta, Rs=Rs, VR=VR, VRe=VRe, f=f, ph=phase(f))
    qs = sorted(set(k[0] for k in recs)); bs = sorted(set(k[1] for k in recs)); ks = sorted(set(k[2] for k in recs))
    print(f"# V(R) points: {len(recs)}  q={qs} beta={bs} kappa={ks}\n")
    print("## per-(q,beta,kappa): phase via V(R) shape  [a=alpha(Coulomb) s=sigma(string) dV=rise]")
    for q in qs:
        print(f"\n q={q}:")
        for b in bs:
            cells = []
            for k in ks:
                r = recs.get((q, b, k))
                if not r: cells.append(f"k{k:g}:----"); continue
                f = r['f']; tag = {0: "Conf", 1: "Coul", 2: "Higg", -1: "fail"}[r['ph']]
                a = f['alpha'] if f else float('nan'); s = f['sigma'] if f else float('nan')
                cells.append(f"k{k:g}:{tag}(a={a:.2f},s={s:.3f})")
            print(f"   b={b}: " + "  ".join(cells))
    # wedge verdict: at the near-wall column, charge-1 UNSCREENED (Coulomb) for q>=5, SCREENED (Higgs) for q<=4
    print("\n## WEDGE VERDICT (V(R), independent of m_gamma) -- charge-1 screening per q across the column:")
    for q in qs:
        col = [recs[(q, b, k)]['ph'] for b in bs for k in ks if (q, b, k) in recs]
        ncoul = sum(p == 1 for p in col); nhig = sum(p == 2 for p in col); ncon = sum(p == 0 for p in col)
        verdict = "COULOMB WEDGE (unscreened)" if ncoul > nhig and ncoul >= max(1, len(col)//2) else \
                  "HIGGS (screened, no wedge)" if nhig >= ncoul else "mixed"
        print(f"   q={q}: Coulomb={ncoul} Higgs={nhig} Confined={ncon} of {len(col)}  -> {verdict}")

    # figure: V(R) small-multiples, one panel per q (curves colored by (beta,kappa))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    nq = len(qs); fig, axes = plt.subplots(1, nq, figsize=(3.2 * nq, 4.0), squeeze=False)
    for ax, q in zip(axes[0], qs):
        for b in bs:
            for k in ks:
                r = recs.get((q, b, k))
                if not r or r['f'] is None: continue
                ph = r['ph']; c = {0: "#c0392b", 1: "#2b6cb0", 2: "#27ae60"}.get(ph, "0.5")
                ax.errorbar(r['Rs'], r['VR'], yerr=r['VRe'], marker="o", ms=3, lw=1.1, color=c, alpha=0.8,
                            label=f"b{b} k{k}")
        ax.set_title(f"q={q}"); ax.set_xlabel("R")
        if q == qs[0]: ax.set_ylabel("V(R)")
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([0],[0], color="#2b6cb0", lw=2, label="Coulomb/wedge (rises, unscreened)"),
                        Line2D([0],[0], color="#27ae60", lw=2, label="Higgs (flat, screened)"),
                        Line2D([0],[0], color="#c0392b", lw=2, label="Confined (linear)")],
               loc="lower center", ncol=3, fontsize=9)
    fig.suptitle("V(R) wedge cross-check: charge-1 static potential per q (rises=Coulomb wedge, flat=Higgs)", fontsize=12)
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    os.makedirs("u1f_campaign_analysis", exist_ok=True)
    p = "u1f_campaign_analysis/vr_wedge_map.png"; fig.savefig(p, dpi=120); print("\nwrote", p)

if __name__ == "__main__":
    main(sys.argv[1:])
