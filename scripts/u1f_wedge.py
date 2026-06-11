#!/usr/bin/env python3
"""Build the photon-mass m^2(beta,kappa) maps per q from the wedge scan, and test for the q>=5 Coulomb wedge
(massless m2<=0 persisting to large kappa)."""
import glob, re, os
import numpy as np
H = re.compile(r"Ls=(\d+)\s+Lt=(\d+)\s+beta=([\d.]+)\s+kappa=([\d.]+)\s+q=(\d+)")
M = re.compile(r"m_gamma=([\-\d.eE]+)\s+m2=([\-\d.eE]+)\s+\+\-\s+([\-\d.eE]+)")
recs = []
for d in ("u1f_pm", "u1f_pm_lenore"):
    for fn in glob.glob(d + "/job_pm_*.out"):
        t = open(fn).read(); h, m = H.search(t), M.search(t)
        if h and m:
            recs.append(dict(beta=float(h[3]), kappa=float(h[4]), q=int(h[5]),
                             mg=float(m[1]), m2=float(m[2]), m2e=float(m[3])))
print(f"# pm points: {len(recs)}")

MASSLESS = 0.005   # m2 <= this (within ~error) => massless (Coulomb/wedge)
print("\n## m2 vs beta at kappa=2.0 (deep Higgs) per q -- the WEDGE test (m2<=0 => massless Coulomb):")
for q in [2, 4, 5, 6, 8]:
    row = sorted([(r['beta'], r['m2'], r['m2e']) for r in recs if r['q'] == q and abs(r['kappa'] - 2.0) < 1e-6])
    if row:
        s = "  ".join(f"b{b:g}:{m2:+.3f}" for b, m2, e in row)
        nmassless = sum(1 for b, m2, e in row if m2 <= MASSLESS)
        print(f"  q={q} k=2.0: {s}   [{nmassless}/{len(row)} massless]")

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
except Exception as e:
    print("no matplotlib:", e); raise SystemExit
os.makedirs("u1f_campaign_analysis", exist_ok=True)
qs = [2, 4, 5, 6, 8]
fig, axes = plt.subplots(1, 5, figsize=(19, 3.8))
for ax, q in zip(axes, qs):
    pts = [r for r in recs if r['q'] == q]
    if not pts: ax.set_title(f"q={q} (no data)"); continue
    bs = sorted(set(r['beta'] for r in pts)); ks = sorted(set(r['kappa'] for r in pts))
    g = np.full((len(ks), len(bs)), np.nan)
    for r in pts: g[ks.index(r['kappa']), bs.index(r['beta'])] = r['m2']
    vmax = np.nanmax(np.abs(g));
    im = ax.imshow(g, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   extent=[min(bs), max(bs), min(ks), max(ks)])
    # m2=0 contour (massless boundary)
    try:
        B, K = np.meshgrid(bs, ks)
        ax.contour(B, K, g, levels=[0.0], colors="k", linewidths=1.5)
    except Exception: pass
    ax.set_xlabel("beta"); ax.set_ylabel("kappa"); ax.set_title(f"q={q}: m^2 (blue=massless=Coulomb)")
    fig.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle("Photon mass m^2(beta,kappa), 16^3x8 -- q>=5 wedge = blue (massless) persisting to large kappa", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
p = "u1f_campaign_analysis/wedge_mphoton.png"; fig.savefig(p, dpi=110); print("\nwrote", p)
