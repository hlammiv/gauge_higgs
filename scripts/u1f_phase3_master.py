#!/usr/bin/env python3
"""MASTER full-plane phase3 (per q) -- stitches two datasets, each authoritative where it is best:
  * Tier A point mesh (L=12, beta 0.85-1.45, kappa 0-0.9): the CONFINED wall + triple point + kappa=0 axis,
    classified by sigma_1 (charge-1 string tension) x <cos>_link (matter). [cheap, dense, reaches kappa=0]
  * Combined pm+V(R) grid (L_s=20, beta 1.1-2.5, kappa>=0.6): the DECONFINED Coulomb/Higgs/wedge, classified
    by sigma_1 x V(R) alpha (the robust same-config discriminant; m_gamma unreliable on the deconfined side).
Overlap (beta 1.1-1.45, kappa 0.6-0.9): the combined layer is drawn on top (it has V(R)). kappa capped at 2.5
for the headline (the wedge kappa->inf asymptote is a separate figure). Different L per region (locating phases,
not scaling -- sigma_1 is L-stable); provenance labeled."""
import glob, re
import numpy as np

H_PT = re.compile(r"point: D=\d+ L=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
H_PM = re.compile(r"pm: Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
COS  = re.compile(r"<cos>_link=([\-\d.eE]+)")
SG1  = re.compile(r"sigma1=([\-\d.eEna]+)")
WL   = re.compile(r"W R=(\d+) T=(\d+)\s+W=([\d.eE+\-]+)\s+\+-\s+([\d.eE+\-]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

import sys; sys.path.insert(0, 'scripts'); import u1f_vr as vr
S1_THR, COS_THR, ALPHA_THR = 0.15, 0.50, 0.03

# --- Tier A (confined/triple-point/kappa=0): sigma_1 x <cos> ---
tierA = {}
for d in ("u1f_fg_lenore", "u1f_fg"):
    for fn in glob.glob(d + "/job_point_*.out"):
        t = open(fn).read(); h = H_PT.search(t)
        if not h or "# DONE" not in t: continue
        s, c = SG1.search(t), COS.search(t)
        tierA[(int(h[4]), float(h[2]), float(h[3]))] = (f(s[1]) if s else np.nan, f(c[1]) if c else np.nan)
# --- Combined (deconfined Coulomb/Higgs): sigma_1 x V(R) alpha ---
comb = {}
alpha = {}   # (q,beta,kappa) -> V(R) Coulomb coeff (combined grid only) -- for the WEDGE overlay
for fn in glob.glob("u1f_fg/job_pm_*.out"):
    t = open(fn).read(); h = H_PM.search(t)
    if not h or "# DONE" not in t: continue
    s, c = SG1.search(t), COS.search(t)
    W = {(int(R), int(T)): (float(w), float(e)) for R, T, w, e in WL.findall(t)}
    Rs, VR, VRe = vr.V_of_R(W, max((R for (R, _T) in W), default=0)); ff = vr.fit(Rs, VR, VRe)
    key = (int(h[5]), float(h[3]), float(h[4]))
    comb[key] = (f(s[1]) if s else np.nan, f(c[1]) if c else np.nan)   # (sigma_1, <cos>) -> same classifier
    alpha[key] = ff['alpha'] if ff else np.nan
# --- L_s=16 base (coarse full-plane backbone; fills the gaps between Tier A and combined): sigma_1 x <cos> ---
base = {}
for d in ("u1f_gauge", "u1f_gauge_lenore"):
    for fn in glob.glob(d + "/*.out"):
        t = open(fn).read(); h = H_PM.search(t)
        if not h: continue
        s, c = SG1.search(t), COS.search(t)
        if not s: continue
        base[(int(h[5]), float(h[3]), float(h[4]))] = (f(s[1]), f(c[1]) if c else np.nan)
print(f"# base(L_s=16) cells: {len(base)}   Tier A cells: {len(tierA)}   combined cells: {len(comb)}")

# CONSISTENT 3-phase classifier EVERYWHERE: sigma_1 (confinement) x <cos> (matter). Monotonic, q-blind matter
# axis -> every q has a Higgs at large kappa. (Earlier V(R)-alpha-primary scheme made the wedge eat the Higgs.)
def cls(s1, cos):
    if np.isfinite(s1) and s1 > S1_THR: return 0          # Confined (incl. deep-Z_q-confined)
    return 2 if (np.isfinite(cos) and cos > COS_THR) else 1  # Higgs (matter condensed) / Coulomb (disordered)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60"])
TRIPLE = {2: (0.95, 0.37), 3: (0.98, 0.50), 4: (1.045, 0.585), 5: (1.047, 0.59), 6: (1.048, 0.59), 8: (1.047, 0.59)}
KMAX = 2.5
qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
def layer(ax, cells, q):
    pts = [(b, k, v) for (qq, b, k), v in cells.items() if qq == q and k <= KMAX + 1e-9]
    if not pts: return
    bs = sorted(set(b for b, k, v in pts)); ks = sorted(set(k for b, k, v in pts))
    PH = np.full((len(ks), len(bs)), np.nan)
    for b, k, v in pts: PH[ks.index(k), bs.index(b)] = cls(*v)
    B, K = np.meshgrid(bs, ks)
    ax.pcolormesh(B, K, np.ma.masked_invalid(PH), cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest")
for ax, q in zip(axes.flat, qs):
    layer(ax, base, q)                # coarse full-plane backbone (fills gaps)
    layer(ax, tierA, q)               # fine confined wall / triple point / kappa=0
    layer(ax, comb, q)                # fine deconfined (same sigma_1 x <cos> classifier -> monotonic, no seam)
    # WEDGE overlay (#26): inside the Higgs (matter condensed) but photon LIGHT (V(R) alpha>0.03) -- q>=5 feature.
    wx = [b for (qq, b, k) in alpha if qq == q and k <= KMAX + 1e-9
          and (qq, b, k) in comb and cls(*comb[(qq, b, k)]) == 2
          and np.isfinite(alpha[(qq, b, k)]) and alpha[(qq, b, k)] > ALPHA_THR]
    wy = [k for (qq, b, k) in alpha if qq == q and k <= KMAX + 1e-9
          and (qq, b, k) in comb and cls(*comb[(qq, b, k)]) == 2
          and np.isfinite(alpha[(qq, b, k)]) and alpha[(qq, b, k)] > ALPHA_THR]
    if wx: ax.scatter(wx, wy, s=70, facecolors="none", edgecolors="cyan", linewidths=1.6, zorder=18,
                      marker="s", label="_wedge")
    if q in TRIPLE:
        ax.plot(*TRIPLE[q], "*", color="gold", ms=20, mec="k", mew=1.0, zorder=20)
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(0.4, 2.5); ax.set_ylim(0, KMAX)
fig.legend(handles=[Patch(facecolor=cmap(i), label=l) for i, l in enumerate(["Confined", "Coulomb", "Higgs"])] +
           [plt.Line2D([0], [0], marker="*", color="gold", mec="k", ms=15, ls="", label="triple point"),
            plt.Line2D([0], [0], marker="s", mfc="none", mec="cyan", mew=1.6, ms=9, ls="",
                       label="Coulomb WEDGE (in Higgs, V(R) photon light; q≥5)")],
           loc="lower center", ncol=5, fontsize=9)
fig.suptitle("U(1)+charge-q Higgs MASTER phase diagram (frozen): Confined / Coulomb / Higgs (σ₁×⟨cos⟩, full plane)\n"
             "every q has a Higgs at large κ; the intermediate-Coulomb WEDGE (q≥5) is the cyan-hatched sub-region inside it",
             fontsize=12)
fig.tight_layout(rect=[0, 0.05, 1, 0.95])
import os; os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/phase3_master.png"; fig.savefig(p, dpi=130); print("wrote", p)
