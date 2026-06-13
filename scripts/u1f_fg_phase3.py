#!/usr/bin/env python3
"""Phase3-style (beta,kappa) per-q map built DIRECTLY from the combined pm+V(R) fine grid (u1f_fg) -- same
config stream, both Coulomb metrics. Classify by the ROBUST discriminants: Confined sigma_1>0.15 (charge-1
area law); else Coulomb if V(R) alpha>0.03 (charge-1 unscreened, the wedge) else Higgs (screened). White =
not-yet-run (kappa refinement still in flight). Preliminary view of the INCOMPLETE grid to judge next steps."""
import sys, glob, re
sys.path.insert(0, 'scripts')
import numpy as np, u1f_vr as vr

H   = re.compile(r"Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
SG1 = re.compile(r"sigma1=([\-\d.eEna]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

recs = {}
for fn in glob.glob("u1f_fg/job_pm_*.out"):
    t = open(fn).read(); h = H.search(t)
    if not h or "# DONE" not in t: continue
    s1 = SG1.search(t)
    meta, W = vr.load(fn); R, VR, VRe = vr.V_of_R(W, meta['Rmax']); ff = vr.fit(R, VR, VRe)
    recs[(int(h[5]), float(h[3]), float(h[4]))] = dict(
        s1=f(s1[1]) if s1 else float('nan'), alpha=ff['alpha'] if ff else float('nan'))
print(f"# combined cells DONE: {len(recs)}")

S1_THR, ALPHA_THR = 0.15, 0.03
def cls(r):
    if np.isfinite(r['s1']) and r['s1'] > S1_THR: return 0          # Confined
    if np.isfinite(r['alpha']) and r['alpha'] > ALPHA_THR: return 1  # Coulomb (V(R) unscreened)
    return 2                                                         # Higgs

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60"])
qs = [2, 3, 4, 5, 6, 8]
betas = sorted(set(k[1] for k in recs)); kaps = sorted(set(k[2] for k in recs))
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ax, q in zip(axes.flat, qs):
    PH = np.full((len(kaps), len(betas)), np.nan)
    for (qq, b, k), r in recs.items():
        if qq != q: continue
        PH[kaps.index(k), betas.index(b)] = cls(r)
    B, K = np.meshgrid(betas, kaps)
    ax.pcolormesh(B, K, np.ma.masked_invalid(PH), cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest")
    ax.set_facecolor("white")
    ncov = np.isfinite(PH).sum()
    ax.set_title(f"q={q}  ({ncov}/{len(betas)*len(kaps)} cells)")
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$")
    ax.set_xticks(betas); ax.set_yticks(kaps); ax.tick_params(labelsize=7)
fig.legend(handles=[Patch(facecolor=cmap(i), label=l) for i, l in enumerate(
           ["Confined (σ₁>0.15)", "Coulomb (V(R) α>0.03)", "Higgs (screened)"])] +
           [Patch(facecolor="white", edgecolor="0.7", label="not yet run")],
           loc="lower center", ncol=4, fontsize=10)
fig.suptitle("Phase3 from the combined pm+V(R) fine grid (u1f_fg): deconfined region, V(R)-α Coulomb/Higgs split.\n"
             "q≥4 full grid (β 1.1–2.5, κ 0.6–7); q2,3 flat Higgs (sampled coarsely — white = unsampled, all Higgs)",
             fontsize=11)
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
import os; os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/phase3_fg_prelim.png"; fig.savefig(p, dpi=120); print("wrote", p)
# coverage report per q
print("\n## coverage (cells done) per q:")
for q in qs:
    ks = sorted(set(k for (qq, b, k) in recs if qq == q))
    bs = sorted(set(b for (qq, b, k) in recs if qq == q))
    print(f"  q={q}: {sum(1 for k in recs if k[0]==q)} cells, beta={bs}, kappa={ks}")
