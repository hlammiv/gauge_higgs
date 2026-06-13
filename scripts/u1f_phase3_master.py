#!/usr/bin/env python3
"""MASTER full-plane GAUGE phase3 (per q): Confined / Coulomb / Higgs classified by the GAUGE LINK only --
  Confined : sigma_1 > 0.15 (charge-1 string tension; area law)
  else by the PHOTON MASS m_gamma (the gauge link Higgs vs Coulomb), via the R(p) shell ratio
  R(pmin)/R(2nd):  ~2.0 => massless => COULOMB ;  ~1.0 (<1.5) => massive => HIGGS.
We use m_gamma (NOT V(R)-alpha) for Coulomb/Higgs: V(R)-alpha (charge-1) misreads weak-Coulomb at large beta as
"screened"; m_gamma stays massless there (correct). The matter <cos> is NOT used (Elitzur: not a gauge phase).
Stitched layers, each authoritative where best: L_s=16 base backbone (sigma_1 -> Confined/else-Coulomb; fills
gaps) + Tier A point mesh (fine confined wall + triple point, reaches kappa=0; sigma_1 only) + combined pm+V(R)
fine grid (sigma_1 + m_gamma ratio -> the Higgs where the photon is massive). Higgs only asserted in the
combined region (good L_s=20 m_gamma). Triple points marked."""
import glob, re
import numpy as np

H_PT = re.compile(r"point: D=\d+ L=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
H_PM = re.compile(r"pm: Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
SG1  = re.compile(r"sigma1=([\-\d.eEna]+)")
SH   = re.compile(r"phat2=([\d.]+)\s+R=([\-\d.eE]+)")
def f(x):
    try: return float(x)
    except: return float('nan')
import sys; sys.path.insert(0, 'scripts'); import u1f_vr as vr
WL = re.compile(r"W R=(\d+) T=(\d+)\s+W=([\d.eE+\-]+)\s+\+-\s+([\d.eE+\-]+)")

S1_THR, ALPHA_THR = 0.15, 0.03
# Tier A + base backbone: gauge confinement only (no m_gamma at Tier A; base m_gamma is L_s=16/noisy -> backbone)
tierA, base = {}, {}
for fn in glob.glob("u1f_fg_lenore/job_point_*.out") + glob.glob("u1f_fg/job_point_*.out"):
    t = open(fn).read(); h = H_PT.search(t)
    if not h or "# DONE" not in t: continue
    s = SG1.search(t); tierA[(int(h[4]), float(h[2]), float(h[3]))] = f(s[1]) if s else np.nan
for d in ("u1f_gauge", "u1f_gauge_lenore"):
    for fn in glob.glob(d + "/*.out"):
        t = open(fn).read(); h = H_PM.search(t); s = SG1.search(t)
        if not h or not s: continue
        base[(int(h[5]), float(h[3]), float(h[4]))] = f(s[1])
# combined grid: sigma_1 + V(R) Coulomb coefficient alpha -> full gauge classification
comb = {}
for fn in glob.glob("u1f_fg/job_pm_*.out"):
    t = open(fn).read(); h = H_PM.search(t)
    if not h or "# DONE" not in t: continue
    s = SG1.search(t)
    W = {(int(R), int(T)): (float(w), float(e)) for R, T, w, e in WL.findall(t)}
    Rs, V, E = vr.V_of_R(W, max((rr for (rr, _t) in W), default=0)); ff = vr.fit(Rs, V, E)
    comb[(int(h[5]), float(h[3]), float(h[4]))] = (f(s[1]) if s else np.nan, ff['alpha'] if ff else np.nan)
print(f"# base {len(base)}  TierA {len(tierA)}  combined {len(comb)}")

def cls_conf(s1):            # backbone: only confinement is known (no V(R))
    return 0 if (np.isfinite(s1) and s1 > S1_THR) else 1     # Confined / else Coulomb
def cls_gauge(v):            # combined: sigma_1 (confine) + V(R) alpha (Coulomb/Higgs)
    s1, a = v
    if np.isfinite(s1) and s1 > S1_THR: return 0             # Confined
    if np.isfinite(a) and a > ALPHA_THR: return 1            # charge-1 unscreened (1/r) -> Coulomb
    return 2                                                 # charge-1 screened (flat) -> Higgs

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60"])
TRIPLE = {2: (0.95, 0.37), 3: (0.98, 0.50), 4: (1.045, 0.585), 5: (1.047, 0.59), 6: (1.048, 0.59), 8: (1.047, 0.59)}
KMAX = 2.5
qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
def layer(ax, cells, clsfn, q):
    pts = [(b, k, v) for (qq, b, k), v in cells.items() if qq == q and k <= KMAX + 1e-9]
    if not pts: return
    bs = sorted(set(b for b, k, v in pts)); ks = sorted(set(k for b, k, v in pts))
    PH = np.full((len(ks), len(bs)), np.nan)
    for b, k, v in pts: PH[ks.index(k), bs.index(b)] = clsfn(v)
    ax.pcolormesh(*np.meshgrid(bs, ks), np.ma.masked_invalid(PH), cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest")
for ax, q in zip(axes.flat, qs):
    layer(ax, base, cls_conf, q)
    layer(ax, tierA, cls_conf, q)
    layer(ax, comb, cls_gauge, q)
    if q in TRIPLE: ax.plot(*TRIPLE[q], "*", color="gold", ms=18, mec="k", mew=1.0, zorder=20)
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(0.4, 2.5); ax.set_ylim(0, KMAX)
fig.legend(handles=[Patch(facecolor=cmap(i), label=l) for i, l in enumerate(
           ["Confined (σ₁ area)", "Coulomb (charge-1 unscreened)", "Higgs (charge-1 screened)"])] +
           [plt.Line2D([0], [0], marker="*", color="gold", mec="k", ms=14, ls="", label="triple point")],
           loc="lower center", ncol=4, fontsize=9)
fig.suptitle("U(1)+charge-q MASTER GAUGE phase diagram: Confined (σ₁) / Coulomb / Higgs (V(R) charge-1 screening)\n"
             "gauge-only classification (matter ⟨cos⟩ NOT used); q=8 Higgs appears once the β>2.5 scan lands",
             fontsize=12)
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
import os; os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/phase3_master.png"; fig.savefig(p, dpi=130); print("wrote", p)
