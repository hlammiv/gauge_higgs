#!/usr/bin/env python3
"""Clean per-q (beta,kappa) 3-phase diagram from the definitive L=8 scan.
Classify with the SHARP order parameters: |P1| (charge-1 Polyakov: ~0 confined, >0 deconfined -- jumps at the
confinement line) and <cos>_link (matter ordering -> Higgs). Draw the two boundary contours + the triple point."""
import glob, re, os
import numpy as np
H = re.compile(r"L=(\d+)\s+beta=([\d.]+)\s+kappa=([\d.]+)\s+q=(\d+)")
COS = re.compile(r"<cos>_link=([\-\d.eE]+)")
P1 = re.compile(r"\|P1\|=([\-\d.eEna]+)")
def f(x):
    try: return float(x)
    except: return float('nan')
recs = []
for d in ("u1f_phase", "u1f_phase_lenore"):
    for fn in glob.glob(d + "/job_point_*.out"):
        t = open(fn).read(); h, c, p = H.search(t), COS.search(t), P1.search(t)
        if h and c and p and int(h[1]) == 8:
            recs.append(dict(b=float(h[2]), k=float(h[3]), q=int(h[4]), cos=f(c[1]), P1=f(p[1])))
print(f"# L=8 points: {len(recs)}")

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
COSc, P1c = 0.5, 0.5          # Higgs / deconfinement thresholds
# phase: 0 Confined, 1 Coulomb, 2 Higgs
def cls(cos, p1):
    if cos > COSc: return 2
    return 1 if p1 > P1c else 0
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60"])   # Confined, Coulomb, Higgs
labels = ["Confined", "Coulomb", "Higgs"]

qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ax, q in zip(axes.flat, qs):
    pts = [r for r in recs if r['q'] == q]
    if len(pts) < 6: ax.set_title(f"q={q} (no data)"); continue
    bs = sorted(set(r['b'] for r in pts)); ks = sorted(set(r['k'] for r in pts))
    COSg = np.full((len(ks), len(bs)), np.nan); P1g = np.full((len(ks), len(bs)), np.nan)
    PH = np.full((len(ks), len(bs)), np.nan)
    for r in pts:
        i, j = ks.index(r['k']), bs.index(r['b'])
        COSg[i, j] = r['cos']; P1g[i, j] = r['P1']; PH[i, j] = cls(r['cos'], r['P1'])
    B, K = np.meshgrid(bs, ks)
    ax.pcolormesh(B, K, PH, cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest", alpha=0.9)
    # confinement line (|P1| jump) -- ONLY in the non-Higgs strip (mask the deep-Zq wiggle inside Higgs)
    P1m = np.where(COSg <= COSc, P1g, np.nan)
    try: ax.contour(B, K, P1m, levels=[P1c], colors="k", linewidths=2.2)
    except Exception: pass
    # Higgs onset (matter crossover) -- the boundary of the green region
    try: ax.contour(B, K, COSg, levels=[COSc], colors="k", linewidths=2.0, linestyles="--")
    except Exception: pass
    # triple point ~ where the confinement line meets the Higgs onset: the beta where |P1|=0.5 at the Higgs-onset kappa
    try:
        # Higgs-onset kappa per beta (interp), and confinement beta per kappa (interp) -> intersection
        kh = []  # (beta, kappa_higgs)
        for j, b in enumerate(bs):
            col = COSg[:, j]
            for i in range(1, len(ks)):
                if col[i-1] < COSc <= col[i]:
                    kh.append((b, ks[i-1] + (COSc-col[i-1])*(ks[i]-ks[i-1])/(col[i]-col[i-1]))); break
        bc = []  # (kappa, beta_conf) from |P1|=0.5 in the low-kappa strip
        for i, k in enumerate(ks):
            row = P1g[i, :]
            for j in range(1, len(bs)):
                if row[j-1] < P1c <= row[j]:
                    bc.append((k, bs[j-1] + (P1c-row[j-1])*(bs[j]-bs[j-1])/(row[j]-row[j-1]))); break
        if kh and bc:
            # find a (beta,kappa) close to both curves
            import numpy as _np
            kharr = _np.array(kh); bcarr = _np.array(bc)
            best = None
            for (b0, k0) in kharr:
                k_at_b = _np.interp(b0, bcarr[:,1][_np.argsort(bcarr[:,1])], bcarr[:,0][_np.argsort(bcarr[:,1])]) if len(bcarr)>1 else None
                if k_at_b is not None and abs(k_at_b - k0) < (best[2] if best else 1e9):
                    best = (b0, k0, abs(k_at_b - k0))
            if best and best[2] < 0.35:
                ax.plot(best[0], best[1], "*", color="gold", ms=20, mec="k", mew=1.0, zorder=20)
    except Exception: pass
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(min(bs), max(bs)); ax.set_ylim(min(ks), 2.0)   # focus 0-2 (Confined/Coulomb live at low kappa)
    if q >= 5:
        ax.text(0.5, 0.92, "(large-$\\kappa$ Higgs hides a Coulomb wedge, m$_\\gamma$=0)",
                transform=ax.transAxes, ha="center", fontsize=7, style="italic", color="0.2")
fig.legend(handles=[Patch(facecolor=cmap(i), label=labels[i]) for i in range(3)] +
           [plt.Line2D([0], [0], color="k", lw=2, label="confinement line (|P$_1$|=0.5)"),
            plt.Line2D([0], [0], color="k", lw=2, ls="--", label="Higgs onset (<cos>=0.5)")],
           loc="lower center", ncol=5, fontsize=9)
fig.suptitle("U(1)+charge-q Higgs phase diagram (frozen, L=8) -- Confined / Coulomb / Higgs", fontsize=14)
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
p = "u1f_campaign_analysis/phase3_clean.png"; fig.savefig(p, dpi=120); print("wrote", p)
