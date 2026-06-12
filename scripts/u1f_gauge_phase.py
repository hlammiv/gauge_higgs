#!/usr/bin/env python3
"""Per-q (beta,kappa) 3-phase diagram classified ENTIRELY by GAUGE-FIELD observables, on the unified
16^3 x 8 grid where EVERY point carries the photon mass m_gamma. This is the figure that puts the massless
photon (Coulomb) region onto the map.

Classification (definitive set {m_gamma, sigma_1, rho_M}; see docs/gauge_phase_observables_DEFINITIVE.md):
  Confined : rho_M high              (monopoles condense; charge-1 AND charge-q area law)
  Coulomb  : not confined AND m_gamma ~ 0   (massless photon -- the UNIQUE Coulomb signature; q>=5 wedge)
  Higgs    : not confined AND m_gamma > 0    (massive photon; matter condensed)
The rho_M and m^2 thresholds are set data-drivenly from the corners (printed below)."""
import glob, re, os
import numpy as np

H   = re.compile(r"beta=([\d.]+)\s+kappa=([\d.]+)\s+q=(\d+)")
M2  = re.compile(r"m2=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
COS = re.compile(r"<cos>_link=([\-\d.eE]+)")
RHO = re.compile(r"rho_M=([\-\d.eE]+)")
P1  = re.compile(r"\|P1\|=([\-\d.eEna]+)")
SG1 = re.compile(r"sigma1=([\-\d.eEna]+)")
SGq = re.compile(r"sigmaq=([\-\d.eEna]+)")
def f(x):
    try: return float(x)
    except: return float('nan')

recs = []
for d in ("u1f_gauge", "u1f_gauge_lenore"):
    for fn in glob.glob(d + "/*.out"):
        t = open(fn).read()
        h, m, rho = H.search(t), M2.search(t), RHO.search(t)
        if not (h and m and rho): continue
        c, p, s1, sq = COS.search(t), P1.search(t), SG1.search(t), SGq.search(t)
        recs.append(dict(b=float(h[1]), k=float(h[2]), q=int(h[3]),
                         m2=f(m[1]), m2e=f(m[2]),
                         rho=f(rho[1]), cos=f(c[1]) if c else float('nan'),
                         P1=f(p[1]) if p else float('nan'),
                         s1=f(s1[1]) if s1 else float('nan'),
                         sq=f(sq[1]) if sq else float('nan')))
# de-dup (lucia+lenore overlap) keep last
seen = {}
for r in recs: seen[(r['q'], r['b'], r['k'])] = r
recs = list(seen.values())
print(f"# gauge points: {len(recs)}   q={sorted(set(r['q'] for r in recs))}")

# ---- data-driven thresholds from the corners ----
def corner(q, b, k):
    for r in recs:
        if r['q']==q and abs(r['b']-b)<1e-6 and abs(r['k']-k)<1e-6: return r
    return None
print("\n## corner survey (rho_M, m^2, sigma1) to set thresholds:")
for (lbl, b, k) in [("CONFINED  beta=0.4 k=0", 0.4, 0.0), ("COULOMB   beta=2.5 k=0", 2.5, 0.0),
                    ("deepHiggs beta=1.0 k=2.5", 1.0, 2.5)]:
    for q in (2, 6):
        r = corner(q, b, k)
        if r: print(f"  q={q} {lbl}:  rho_M={r['rho']:8.2f}  m2={r['m2']:+.4f}+-{r['m2e']:.4f}  s1={r['s1']:.3f}  |P1|={r['P1']:.3f}")

# TWO-AXIS classification (the validated SU(2) recipe) on the ROBUST channels:
#   GAUGE axis  = sigma_1 (charge-1 Wilson string tension). SHARP. Locates confinement: at kappa=0 it drops
#                 across beta~1.0 = the pure-gauge confinement->Coulomb transition (beta_c~1.01). A genuine
#                 gauge observable. (sigma_q, the charge-q Creutz ratio, is NOISE at these loop sizes -- nan or
#                 garbage-NEGATIVE even in the deconfined Coulomb phase -- so it was wrong to use it; doing so
#                 ate the entire kappa=0 Coulomb phase. sigma_1 is clean.)
#   MATTER axis = <cos>_link. CROSSOVER (Elitzur: the Higgs has NO local gauge order parameter). =0 at kappa=0
#                 (matter decoupled), rises monotonically; q-BLIND (condenses at the same kappa for all q).
# Phases:
#   Higgs    : matter condensed (<cos> high) -- regardless of whether residual Z_q still confines charge-1.
#   Confined : matter disordered AND charge-1 confined (sigma_1 high).
#   Coulomb  : matter disordered AND charge-1 free (sigma_1 low). At kappa=0 this is the large-beta photon.
# NOTE: this 2-axis map CANNOT see the q>=5 Coulomb WEDGE -- the wedge is Coulomb persisting INSIDE the
# matter-condensed region, and <cos> is q-blind (condenses identically for all q). Only m_gamma can see it,
# and m_gamma is unreliable at L_s=16 (degenerate R-fits). We OVERLAY m_gamma-massless cells as tentative
# wedge candidates, but the wedge needs the dedicated L_s>=20 m_gamma run.
S1_THR  = 0.15        # charge-1 area-law cut (deconf sigma_1 ~ 0-0.09; conf ~ 0.35-1.6) -- sharp
COS_THR = 0.50        # matter-condensation (Higgs) crossover; soft (Elitzur), threshold-dependent
M2_THR  = 0.013       # massless cut for the tentative wedge overlay (~2x median R-fit error)
print(f"\n# thresholds: sigma_1> {S1_THR} charge-1 CONFINED (sharp) ; <cos>> {COS_THR} matter CONDENSED (crossover)")

# 0 Confined, 1 Coulomb, 2 Higgs
def cls(r):
    matter = np.isfinite(r['cos']) and r['cos'] > COS_THR    # Higgs: matter condensed
    s1conf = np.isfinite(r['s1'])  and r['s1']  > S1_THR     # charge-1 confined (sharp gauge axis)
    if matter:   return 2     # Higgs (matter condensed)
    if s1conf:   return 0     # Confined (charge-1 confined, matter disordered)
    return 1                  # Coulomb (charge-1 free, matter disordered)

def wedge_candidate(r):
    # inside the Higgs, charge-1 free, photon CONFIDENTLY massless -> tentative Coulomb-wedge cell (m_gamma noisy!)
    return (cls(r) == 2 and np.isfinite(r['s1']) and r['s1'] < S1_THR
            and np.isfinite(r['m2']) and 0.0 <= r['m2'] <= M2_THR)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60"])   # Confined, Coulomb, Higgs
labels = ["Confined", "Coulomb (m$_\\gamma$=0)", "Higgs (matter condensed)"]

qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ax, q in zip(axes.flat, qs):
    pts = [r for r in recs if r['q'] == q]
    if len(pts) < 6: ax.set_title(f"q={q} (no data)"); continue
    bs = sorted(set(r['b'] for r in pts)); ks = sorted(set(r['k'] for r in pts))
    PH = np.full((len(ks), len(bs)), np.nan)
    S1g = np.full((len(ks), len(bs)), np.nan); COSg = np.full((len(ks), len(bs)), np.nan)
    for r in pts:
        i, j = ks.index(r['k']), bs.index(r['b'])
        PH[i, j] = cls(r); S1g[i, j] = r['s1']; COSg[i, j] = r['cos']
    B, K = np.meshgrid(bs, ks)
    ax.pcolormesh(B, K, PH, cmap=cmap, vmin=-0.5, vmax=2.5, shading="nearest", alpha=0.85)
    # SHARP confinement line: sigma_1 = S1_THR (the genuine gauge observable -- locates beta_c at kappa=0)
    try: ax.contour(B, K, S1g, levels=[S1_THR], colors="k", linewidths=2.2)
    except Exception: pass
    # SOFT Higgs onset: <cos> = COS_THR (Elitzur crossover -- dashed to flag it is not a sharp transition)
    try: ax.contour(B, K, COSg, levels=[COS_THR], colors="k", linewidths=1.8, linestyles="--")
    except Exception: pass
    # tentative wedge candidates (m_gamma-massless inside Higgs) -- m_gamma noisy at L_s=16
    wx = [r['b'] for r in pts if wedge_candidate(r)]; wy = [r['k'] for r in pts if wedge_candidate(r)]
    if wx: ax.plot(wx, wy, "o", mfc="none", mec="cyan", mew=1.8, ms=11, zorder=20)
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(min(bs), max(bs)); ax.set_ylim(min(ks), max(ks))
fig.legend(handles=[Patch(facecolor=cmap(i), label=labels[i]) for i in range(3)] +
           [plt.Line2D([0],[0], color="k", lw=2.2, label=r"confinement line $\sigma_1$=%.2f (sharp)" % S1_THR),
            plt.Line2D([0],[0], color="k", lw=1.8, ls="--", label=r"Higgs onset $\langle\cos\rangle$=%.1f (crossover)" % COS_THR),
            plt.Line2D([0],[0], marker="o", mfc="none", mec="cyan", mew=1.8, ms=10, ls="",
                       label=r"wedge candidate (m$_\gamma\!\approx$0 in Higgs; tentative @L$_s$=16)")],
           loc="lower center", ncol=3, fontsize=9)
fig.suptitle("U(1)+charge-q Higgs (frozen, 16$^3\\times$8): $\\sigma_1$ confinement axis (sharp) $\\times$ $\\langle\\cos\\rangle$ matter axis (crossover)",
             fontsize=13)
fig.tight_layout(rect=[0, 0.06, 1, 0.96])
os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/phase3_gauge.png"; fig.savefig(p, dpi=120); print("\nwrote", p)

# ---- per-q phase census + kappa=0 sanity row + wedge-candidate count ----
print("\n## phase census (kappa=0 must show BOTH Confined and Coulomb):")
for q in qs:
    pts = [r for r in recs if r['q'] == q]
    if not pts: continue
    n = [0,0,0]
    for r in pts: n[cls(r)] += 1
    row0 = sorted([r for r in pts if abs(r['k']) < 1e-6], key=lambda r: r['b'])
    seq = "".join("CcH"[cls(r)] for r in row0)   # kappa=0 row, low->high beta
    nw = sum(1 for r in pts if wedge_candidate(r))
    print(f"  q={q}: Confined={n[0]:2d} Coulomb={n[1]:2d} Higgs={n[2]:2d} | kappa=0 row(lo->hi beta)='{seq}' | wedge-cand={nw}")
