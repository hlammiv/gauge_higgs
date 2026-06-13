#!/usr/bin/env python3
"""Does the Coulomb-Higgs wedge edge asymptote as kappa->inf? Plot beta_CH(kappa) per q = the largest beta
still Coulomb (V(R) alpha>0.03) at each kappa, from the combined pm+V(R) grid (now kappa to 7). Flat at high
kappa = asymptoted. usage: scripts/u1f_wedge_asymptote.py"""
import sys, glob
sys.path.insert(0, 'scripts')
import numpy as np, u1f_combined as C
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

rows = C.load_rows(glob.glob("u1f_fg/job_pm_*.out"))
QS = [4, 5, 6, 8]; ALPHA_THR = 0.03
fig, ax = plt.subplots(figsize=(8, 5.4))
for q in QS:
    kaps = sorted(set(k for (qq, b, k) in rows if qq == q))
    bk = []
    for k in kaps:
        bb = [b for (qq, b, kk) in rows if qq == q and kk == k
              and np.isfinite(rows[(q, b, k)]['alpha']) and rows[(q, b, k)]['alpha'] > ALPHA_THR]
        if bb: bk.append((k, max(bb)))
    if bk:
        ks, bs = zip(*bk); ax.plot(ks, bs, "o-", lw=2, ms=7, label=f"q={q}")
ax.set_xlabel(r"$\kappa$"); ax.set_ylabel(r"Coulomb-Higgs boundary $\beta_{CH}$ (wedge edge, V(R) $\alpha$>0.03)")
ax.set_title("Wedge edge asymptote: $\\beta_{CH}(\\kappa)$ per q (κ to 7)\nflat at high κ = asymptoted to κ→∞ value")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout(); p = "u1f_campaign_analysis/wedge_asymptote.png"; fig.savefig(p, dpi=120); print("wrote", p)
for q in QS:
    kaps = sorted(set(k for (qq, b, k) in rows if qq == q))
    print(f"q={q}: " + " ".join(f"k{k:g}:{max([b for (qq,b,kk) in rows if qq==q and kk==k and rows[(q,b,k)]['alpha']>ALPHA_THR], default=0):.2f}" for k in kaps))
