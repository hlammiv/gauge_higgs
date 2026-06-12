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

# Textbook FS classification by Wilson-loop CHARGE CONTENT + photon mass (rho_M NOT used -- it leaks into
# the deep-Higgs-Z_q corner because charge-1 is confined there by the residual Z_q).
#   sigma_q (charge-q, the dynamical matter charge): area-law -> charge-q CONFINED.  Creutz chi is nan or
#       negative or large when the loop is deeply confined; small & non-negative when screened.
#   sigma_1 (charge-1): area-law -> charge-1 CONFINED (genuine Confined, OR residual-Z_q inside the Higgs).
#   m_gamma=0 (m2<=thr): massless photon -- only in a genuine Coulomb phase (all charges free + gapless).
M2_THR = 0.013        # massless cut (~2x median R-fit error)
S1_THR = 0.15         # charge-1 area-law cut (deconf sigma_1 ~ 0-0.07; conf ~ 0.3-1.6)
SQ_THR = 0.20         # charge-q area-law cut (screened ~ 0-0.13; conf = nan/neg/large)
print(f"\n# thresholds: m2<= {M2_THR} massless ; sigma_1> {S1_THR} charge-1 conf ; |sigma_q|> {SQ_THR} (or nan) charge-q conf")

# 0 Confined, 1 Coulomb, 2 Higgs, 3 deconfined-but-m_gamma-UNRESOLVED (honest: L_s=16 floor)
# m_gamma at L_s=16 is unreliable (degenerate R-fits -> negative m^2, false +-0.000, +-88 blowups). So in the
# both-charges-screened region we ONLY split Coulomb/Higgs where m^2 is CLEAN: a confident massless (m2<=thr
# AND not an outlier) -> Coulomb; a confident, significant, sane massive (m2>MASS_HI, |m2|<2, m2>3*err) -> Higgs;
# everything else -> "unresolved" (the wedge needs L_s>=20 + better stats, not a guess).
MASS_HI = 0.06
def cls(r):
    sq, s1, m2, e = r['sq'], r['s1'], r['m2'], r['m2e']
    sq_conf = (not np.isfinite(sq)) or abs(sq) > SQ_THR      # charge-q area law -> Confined
    s1_conf = (np.isfinite(s1) and s1 > S1_THR)             # charge-1 area law (residual Z_q inside Higgs)
    if sq_conf:                       return 0   # charge-q confined  -> Confined (strong coupling)
    if s1_conf:                       return 2   # charge-q screened, charge-1 Z_q-confined -> deep-Higgs (Higgs)
    # both charges screened: clean Coulomb/Higgs only where m_gamma is trustworthy
    if not np.isfinite(m2) or abs(m2) > 2.0:        return 3   # fit blowup -> unresolved
    if 0.0 <= m2 <= M2_THR:                          return 1   # clean massless -> Coulomb
    if m2 > MASS_HI and m2 > 3.0 * max(e, 1e-6):     return 2   # clean, significant massive -> Higgs
    return 3                                                    # marginal / negative-noise -> unresolved

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
cmap = ListedColormap(["#c0392b", "#2b6cb0", "#27ae60", "#bfbfbf"])   # Confined, Coulomb, Higgs, unresolved
labels = ["Confined (charge-q area law)", "Coulomb (m$_\\gamma$=0)", "Higgs (charge-q screened)",
          "deconf., m$_\\gamma$ unresolved @ L$_s$=16"]

qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
for ax, q in zip(axes.flat, qs):
    pts = [r for r in recs if r['q'] == q]
    if len(pts) < 6: ax.set_title(f"q={q} (no data)"); continue
    bs = sorted(set(r['b'] for r in pts)); ks = sorted(set(r['k'] for r in pts))
    PH = np.full((len(ks), len(bs)), np.nan)
    M2g = np.full((len(ks), len(bs)), np.nan)
    for r in pts:
        i, j = ks.index(r['k']), bs.index(r['b'])
        PH[i, j] = cls(r); M2g[i, j] = r['m2']
    B, K = np.meshgrid(bs, ks)
    ax.pcolormesh(B, K, PH, cmap=cmap, vmin=-0.5, vmax=3.5, shading="nearest", alpha=0.9)
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(min(bs), max(bs)); ax.set_ylim(min(ks), max(ks))
fig.legend(handles=[Patch(facecolor=cmap(i), label=labels[i]) for i in range(4)],
           loc="lower center", ncol=4, fontsize=9)
fig.suptitle("U(1)+charge-q Higgs phase diagram (frozen, GAUGE-classified, 16$^3\\times$8) -- m$_\\gamma$ puts the Coulomb region on the map",
             fontsize=13)
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/phase3_gauge.png"; fig.savefig(p, dpi=120); print("\nwrote", p)

# ---- per-q phase census + Coulomb-at-large-kappa (the wedge) ----
print("\n## phase census + massless fraction at kappa>=1.5 (the q>=5 Coulomb wedge):")
for q in qs:
    pts = [r for r in recs if r['q'] == q]
    if not pts: continue
    n = [0,0,0,0]
    for r in pts: n[cls(r)] += 1
    print(f"  q={q}: Confined={n[0]:2d} Coulomb={n[1]:2d} Higgs={n[2]:2d} Unresolved={n[3]:2d}")
