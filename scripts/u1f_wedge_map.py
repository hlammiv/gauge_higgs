#!/usr/bin/env python3
"""beta-resolved wedge map: V(R) Coulomb coefficient alpha(q,beta) (kappa-averaged) over the combined
pm+V(R) fine grid. alpha>0.03 = Coulomb wedge (charge-1 unscreened), ~0 = Higgs (screened). Shows the wedge
as a wall-hugging staircase that closes from low q up as beta increases. Reuses scripts/u1f_combined.load_rows.
usage: scripts/u1f_wedge_map.py [job_pm_*.out ...]"""
import sys, glob
sys.path.insert(0, 'scripts')
import numpy as np, u1f_combined as C

rows = C.load_rows(sys.argv[1:] or glob.glob("u1f_fg/job_pm_*.out"))
qs = [2, 3, 4, 5, 6, 8]; betas = sorted(set(k[1] for k in rows)); KAPS = (1.0, 1.5, 2.5)
A = np.full((len(qs), len(betas)), np.nan)
for i, q in enumerate(qs):
    for j, b in enumerate(betas):
        v = [rows[(q, b, k)]['alpha'] for k in KAPS if (q, b, k) in rows]
        if v: A[i, j] = np.mean(v)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
fig, ax = plt.subplots(figsize=(max(8, 1.3 * len(betas) + 4), 5.5))
im = ax.imshow(A, origin="lower", aspect="auto", cmap="RdYlBu",
               norm=TwoSlopeNorm(vmin=0, vcenter=0.03, vmax=0.10),
               extent=[-0.5, len(betas) - 0.5, -0.5, len(qs) - 0.5])
ax.set_yticks(range(len(qs))); ax.set_yticklabels(qs); ax.set_ylabel("q")
ax.set_xticks(range(len(betas))); ax.set_xticklabels([f"{b:g}" for b in betas]); ax.set_xlabel(r"$\beta$")
for i in range(len(qs)):
    for j in range(len(betas)):
        if np.isfinite(A[i, j]):
            ax.text(j, i, f"{A[i, j]:.02f}", ha="center", va="center", fontsize=8,
                    color="white" if A[i, j] > 0.05 else "black")
fig.colorbar(im, ax=ax, label=r"V(R) Coulomb coeff $\alpha$ (κ-avg)")
ax.set_title("β-resolved wedge map α(q,β): Blue=Coulomb wedge (α>0.03), Red=Higgs\n"
             "wedge is a wall-hugging staircase — closes from low q up as β grows")
fig.tight_layout(); p = "u1f_campaign_analysis/wedge_qbeta_map.png"; fig.savefig(p, dpi=120); print("wrote", p)
print("\n## wedge open (alpha>0.03) beta-range per q:")
for i, q in enumerate(qs):
    ob = [betas[j] for j in range(len(betas)) if np.isfinite(A[i, j]) and A[i, j] > 0.03]
    print(f"  q={q}: {ob if ob else 'never (always Higgs)'}")
