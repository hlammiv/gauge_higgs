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

# --- FOLD IN the reliable L_s=20 m_gamma verdict on the deconfined side (the RESOLVED Coulomb wedge) ---
# On the deconfined (charge-1 free) side, the photon mass is the TRUE Coulomb/Higgs discriminant (<cos> is a
# q-blind crossover that cannot see the wedge). The L_s=20 run (scripts/u1f_wedge20.py) resolved m_gamma there
# via the R(p)-shell ratio: massless->Coulomb (the wedge), massive->Higgs. Prefer it wherever available.
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u1f_wedge20 as w20
W20 = {}    # (q, beta, kappa) -> "massless" | "massive" | "amb"
for (q, b, k), rr in w20.load(("u1f_wedge20", "u1f_wedge20_lenore")).items():
    W20[(q, round(b, 2), round(k, 2))] = w20.state(rr)[0]
print(f"# folded L_s=20 m_gamma verdicts: {len(W20)} points "
      f"(massless={sum(v=='massless' for v in W20.values())} massive={sum(v=='massive' for v in W20.values())})")

# --- HARDENED: prefer the L_s=24 m_gamma verdict (finer floor 0.068, finite-T checked) over L_s=20 ---
H24 = re.compile(r"Ls=(\d+) Lt=(\d+) beta=([\d.]+) kappa=([\d.]+) q=(\d+)")
M24 = re.compile(r"m2=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
S24 = re.compile(r"phat2=([\d.]+)\s+R=([\-\d.eE]+)\s+\+-\s+([\-\d.eE]+)")
def _ratio_state(m2, e, sh):    # R(pmin)/R(2nd): ~2.0 massless, ~1.0 massive; m2 significance breaks ties
    if len(sh) < 2 or not np.isfinite(m2) or not sh[1][1]: return "amb"
    rr = sh[0][1] / sh[1][1]
    if rr > 1.7 and abs(m2) < 0.02: return "massless"
    if rr < 1.45 and m2 > 0.03 and e and m2 > 3 * e: return "massive"
    return "amb"
W24 = {}    # (q,beta,kappa) -> verdict, L_s=24 L_t=8 main grid
for fn in glob.glob("u1f_wedge24/*.out"):
    t = open(fn).read(); h = H24.search(t)
    if not h or int(h[2]) != 8: continue
    sh = sorted([(f(a), f(b), f(c)) for a, b, c in S24.findall(t) if f(b) > 0])
    m = M24.search(t)
    W24[(int(h[5]), round(float(h[3]), 2), round(float(h[4]), 2))] = _ratio_state(f(m[1]), f(m[2]), sh)
print(f"# L_s=24 m_gamma verdicts: {len(W24)} (massless={sum(v=='massless' for v in W24.values())} "
      f"massive={sum(v=='massive' for v in W24.values())})")

# --- independent V(R) static-potential verdict (orthogonal observable): screened(Higgs) vs unscreened(Coulomb) ---
VRV = {}
try:
    import u1f_vr as _vr
    for fn in glob.glob("u1f_vr/job_*.out"):
        meta, W = _vr.load(fn); Rs, VR, VRe = _vr.V_of_R(W, meta['Rmax']); ff = _vr.fit(Rs, VR, VRe)
        if ff is None: continue
        scr = (ff['sigma'] <= _vr.SIG_THR) and (ff['alpha'] <= _vr.ALPHA_THR)   # flat V -> charge-1 screened
        VRV[(meta['q'], round(meta['beta'], 2), round(meta['kappa'], 2))] = "screened" if scr else "unscreened"
    print(f"# V(R) verdicts: {len(VRV)} (unscreened={sum(v=='unscreened' for v in VRV.values())} "
          f"screened={sum(v=='screened' for v in VRV.values())})")
except Exception as e:
    print("# V(R) load skipped:", e)

# NEAR-WALL only: m_gamma discriminates Coulomb/Higgs only where a Higgs mass would be RESOLVABLE (beta<=BETA_MG);
# at weak coupling the Higgs mass m^2*a^2 sits below the floor and reads false-massless for every q.
BETA_MG = 1.6
def mg_verdict(r):
    if r['b'] > BETA_MG: return None
    k = (r['q'], round(r['b'], 2), round(r['k'], 2))
    v = W24.get(k)                                            # prefer hardened L_s=24
    if v in ("massless", "massive"): return v
    v = W20.get(k)                                            # fall back to L_s=20
    return v if v in ("massless", "massive") else None

# 0 Confined, 1 Coulomb, 2 Higgs. CONFINEMENT (sigma_1) FIRST; deconfined side split by the RESOLVED near-wall
# m_gamma (massless=Coulomb wedge, massive=Higgs); else <cos> matter crossover.
def cls(r):
    if np.isfinite(r['s1']) and r['s1'] > S1_THR: return 0    # Confined (connects to kappa=inf)
    v = mg_verdict(r)
    if v == "massless": return 1                             # Coulomb -- RESOLVED wedge
    if v == "massive":  return 2                             # Higgs -- RESOLVED
    return 2 if (np.isfinite(r['cos']) and r['cos'] > COS_THR) else 1   # fallback <cos>

def has_mg(r):  return mg_verdict(r) in ("massless", "massive")
def vr_agree(r):
    """Does the INDEPENDENT V(R) verdict match the phase label? True/False/None(no data or confined)."""
    vv = VRV.get((r['q'], round(r['b'], 2), round(r['k'], 2)))
    if vv is None: return None
    ph = cls(r)
    if ph == 1: return vv == "unscreened"    # Coulomb wedge <-> V(R) charge-1 unscreened
    if ph == 2: return vv == "screened"      # Higgs        <-> V(R) charge-1 screened
    return None

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
    # SOFT Higgs onset: <cos> = COS_THR, ONLY on the deconfined side (it is a Coulomb|Higgs boundary there;
    # inside the confined region matter also condenses but it is NOT a transition -- FS analytic connection).
    COSm = np.where(S1g <= S1_THR, COSg, np.nan)
    try: ax.contour(B, K, COSm, levels=[COS_THR], colors="k", linewidths=1.8, linestyles="--")
    except Exception: pass
    # dots: Coulomb/Higgs label RESOLVED by m_gamma (L_s=24 preferred); markers: INDEPENDENT V(R) check
    dx = [r['b'] for r in pts if has_mg(r)]; dy = [r['k'] for r in pts if has_mg(r)]
    if dx: ax.plot(dx, dy, ".", color="k", ms=4, alpha=0.5, zorder=20)
    ax2x = [r['b'] for r in pts if vr_agree(r) is True];  ax2y = [r['k'] for r in pts if vr_agree(r) is True]
    dgx = [r['b'] for r in pts if vr_agree(r) is False];  dgy = [r['k'] for r in pts if vr_agree(r) is False]
    if ax2x: ax.plot(ax2x, ax2y, "P", mfc="w", mec="k", mew=1.3, ms=9, zorder=21)   # V(R) AGREES
    if dgx:  ax.plot(dgx, dgy, "X", mfc="magenta", mec="k", mew=1.0, ms=10, zorder=21)  # V(R) DISAGREES
    ax.set_xlabel(r"$\beta$"); ax.set_ylabel(r"$\kappa$"); ax.set_title(f"q = {q}")
    ax.set_xlim(min(bs), max(bs)); ax.set_ylim(min(ks), max(ks))
fig.legend(handles=[Patch(facecolor=cmap(i), label=labels[i]) for i in range(3)] +
           [plt.Line2D([0],[0], color="k", lw=2.2, label=r"confinement line $\sigma_1$=0.15 (sharp)"),
            plt.Line2D([0],[0], color="k", lw=1.8, ls="--", label=r"Higgs onset $\langle\cos\rangle$=0.5 (crossover)"),
            plt.Line2D([0],[0], marker=".", color="k", ms=8, ls="", alpha=0.5, label=r"label set by m$_\gamma$ (L$_s$=24 near-wall)"),
            plt.Line2D([0],[0], marker="P", mfc="w", mec="k", mew=1.3, ms=9, ls="", label="V(R) AGREES (independent)"),
            plt.Line2D([0],[0], marker="X", mfc="magenta", mec="k", ms=9, ls="", label="V(R) disagrees")],
           loc="lower center", ncol=3, fontsize=9)
fig.suptitle("U(1)+charge-q Higgs (frozen): Confined ($\\sigma_1$) | Coulomb | Higgs -- deconfined split by HARDENED L$_s$=24 m$_\\gamma$, cross-checked by V(R)",
             fontsize=12)
fig.tight_layout(rect=[0, 0.06, 1, 0.96])
os.makedirs("u1f_campaign_analysis", exist_ok=True)
p = "u1f_campaign_analysis/phase3_gauge.png"; fig.savefig(p, dpi=120); print("\nwrote", p)

# ---- per-q phase census + kappa=0 / kappa=max confined-boundary sanity ----
print("\n## census + topology (kappa=0 row: Confined->Coulomb ; kappa=max: Confined must persist = connects to inf):")
ZQC = {2:0.44, 3:0.70, 4:0.95, 5:1.02, 6:1.05, 8:1.08}   # approx pure-Z_q beta_c (from matching)
for q in qs:
    pts = [r for r in recs if r['q'] == q]
    if not pts: continue
    n = [0,0,0]
    for r in pts: n[cls(r)] += 1
    row0 = sorted([r for r in pts if abs(r['k']) < 1e-6], key=lambda r: r['b'])
    seq0 = "".join("CcH"[cls(r)] for r in row0)
    kmax = max(r['k'] for r in pts)
    rowT = sorted([r for r in pts if abs(r['k']-kmax) < 1e-6], key=lambda r: r['b'])
    confT = [r['b'] for r in rowT if cls(r) == 0]
    bbound = max(confT) if confT else None
    print(f"  q={q}: Conf={n[0]:2d} Coul={n[1]:2d} Higgs={n[2]:2d} | k=0:'{seq0}' | "
          f"k={kmax:g} confined up to beta={bbound} (pure-Z{q} beta_c~{ZQC.get(q,'?')})")
