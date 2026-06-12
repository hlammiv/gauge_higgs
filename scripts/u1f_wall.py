#!/usr/bin/env python3
"""Fine-grained confinement-wall + triple-point map from the Tier-A `point` mesh (L=12, dense beta,kappa).
Per q: phase by sigma_1 (confinement, sharp) x <cos>_link (matter crossover); overlay the sigma_1=0.15 wall,
the <cos>=0.5 Higgs onset, the chi_plaq peak ridge, and mark the TRIPLE POINT (wall meets onset).
usage: scripts/u1f_wall.py [dir ...]   (default u1f_fg_lenore u1f_fg)"""
import sys, glob, re
import numpy as np

H   = re.compile(r"point: D=\d+ L=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
COS = re.compile(r"<cos>_link=([\-\d.eE]+)")
SG1 = re.compile(r"sigma1=([\-\d.eEna]+)")
RHO = re.compile(r"rho_M=([\-\d.eE]+)")
CHP = re.compile(r"chi_plaq=([\-\d.eE]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

dirs = sys.argv[1:] or ["u1f_fg_lenore", "u1f_fg"]
recs = {}
for d in dirs:
    for fn in glob.glob(d + "/job_point_*.out"):
        t = open(fn).read(); h = H.search(t)
        if not h or "# DONE" not in t: continue
        c, s, r, cp = COS.search(t), SG1.search(t), RHO.search(t), CHP.search(t)
        recs[(int(h[4]), float(h[2]), float(h[3]))] = dict(
            cos=f(c[1]) if c else float('nan'), s1=f(s[1]) if s else float('nan'),
            rho=f(r[1]) if r else float('nan'), chip=f(cp[1]) if cp else float('nan'))
print(f"# fine-mesh point records: {len(recs)}  q={sorted(set(k[0] for k in recs))}")

S1_THR, COS_THR = 0.15, 0.50
def cls(r):
    if np.isfinite(r['s1']) and r['s1'] > S1_THR: return 0          # Confined
    return 2 if (np.isfinite(r['cos']) and r['cos'] > COS_THR) else 1  # Higgs / Coulomb

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60"])
qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
triples = {}
for ax, q in zip(axes.flat, qs):
    pts = [(b, k, r) for (qq, b, k), r in recs.items() if qq == q]
    if len(pts) < 6: ax.set_title(f"q={q} (no data)"); continue
    bs = sorted(set(b for b, k, r in pts)); ks = sorted(set(k for b, k, r in pts))
    PH = np.full((len(ks), len(bs)), np.nan); S1g = np.full_like(PH, np.nan)
    COSg = np.full_like(PH, np.nan); CHg = np.full_like(PH, np.nan)
    for b, k, r in pts:
        i, j = ks.index(k), bs.index(b)
        PH[i, j] = cls(r); S1g[i, j] = r['s1']; COSg[i, j] = r['cos']; CHg[i, j] = r['chip']
    B, K = np.meshgrid(bs, ks)
    ax.pcolormesh(B, K, PH, cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest", alpha=0.8)
    cw = co = None
    try: cw = ax.contour(B, K, S1g, levels=[S1_THR], colors="k", linewidths=2.4)         # confinement wall
    except Exception: pass
    try: co = ax.contour(B, K, COSg, levels=[COS_THR], colors="k", linewidths=1.8, linestyles="--")  # Higgs onset
    except Exception: pass
    try: ax.contour(B, K, CHg, levels=[np.nanmax(CHg) * 0.6], colors="orange", linewidths=1.2, alpha=0.7)  # chi_plaq ridge
    except Exception: pass
    # triple point: intersection of the two contours (nearest pair of contour vertices)
    try:
        wp = np.vstack([seg for c in cw.allsegs for seg in c]) if cw and cw.allsegs else None
        op = np.vstack([seg for c in co.allsegs for seg in c]) if co and co.allsegs else None
        if wp is not None and op is not None and len(wp) and len(op):
            d = np.linalg.norm(wp[:, None, :] - op[None, :, :], axis=2)
            i, j = np.unravel_index(np.argmin(d), d.shape)
            tp = 0.5 * (wp[i] + op[j])
            if d[i, j] < 0.12:
                ax.plot(*tp, "*", color="gold", ms=22, mec="k", mew=1.2, zorder=20)
                triples[q] = (round(tp[0], 3), round(tp[1], 3))
    except Exception: pass
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(min(bs), max(bs)); ax.set_ylim(min(ks), max(ks))
fig.legend(handles=[Patch(facecolor=cmap(i), label=l) for i, l in enumerate(["Confined", "Coulomb", "Higgs"])] +
           [plt.Line2D([0],[0], color="k", lw=2.4, label=r"confinement wall $\sigma_1$=0.15"),
            plt.Line2D([0],[0], color="k", lw=1.8, ls="--", label=r"Higgs onset $\langle\cos\rangle$=0.5"),
            plt.Line2D([0],[0], color="orange", lw=1.2, label=r"$\chi_{plaq}$ ridge"),
            plt.Line2D([0],[0], marker="*", color="gold", mec="k", ms=14, ls="", label="triple point")],
           loc="lower center", ncol=4, fontsize=9)
fig.suptitle("Fine-grained confinement wall + triple point (Tier A, point L=12, $\\Delta\\beta$=0.05 $\\Delta\\kappa$=0.1)", fontsize=13)
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
import os; os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/wall_triplepoint.png"; fig.savefig(p, dpi=120); print("wrote", p)
print("\n## triple-point locations (beta, kappa) per q:")
for q in qs:
    print(f"  q={q}: {triples.get(q, 'not bracketed in this window')}")
