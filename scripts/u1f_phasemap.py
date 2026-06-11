#!/usr/bin/env python3
"""Classify each (beta,kappa) point of the definitive frozen-U(1) scan into Confined/Coulomb/Higgs using the
DEFINITIVE gauge discriminants (docs/gauge_phase_observables_DEFINITIVE.md), and draw per-q phase maps."""
import glob, re, os
import numpy as np
H = re.compile(r"L=(\d+)\s+beta=([\d.]+)\s+kappa=([\d.]+)\s+q=(\d+)")
OBS = re.compile(r"<cos>_link=([\-\d.eE]+)")
DISC = re.compile(r"rho_M=([\-\d.eEna]+)\s+\|P1\|=([\-\d.eEna]+)\s+\|Pq\|=([\-\d.eEna]+)\s+sigma1=([\-\d.eEna]+)\s+sigmaq=([\-\d.eEna]+)")

def f(x):
    try: return float(x)
    except: return float('nan')

recs = []
for d in ("u1f_phase", "u1f_phase_lenore"):
    for fn in glob.glob(d + "/job_point_*.out"):
        t = open(fn).read(); h, o, dc = H.search(t), OBS.search(t), DISC.search(t)
        if h and o and dc:
            recs.append(dict(L=int(h[1]), beta=float(h[2]), kappa=float(h[3]), q=int(h[4]),
                             cos=f(o[1]), rho=f(dc[1]), P1=f(dc[2]), Pq=f(dc[3]), sig1=f(dc[4]), sigq=f(dc[5])))
print(f"# classified points: {len(recs)}")

# thresholds
COS_ORD = 0.45    # matter ordered (Higgs) -- crossover, soft
SIG_CONF = 0.20   # charge-1 area law (confined)
SIGQ_SCR = 0.05   # charge-q screened (~perimeter)

# phase codes: 0 Coulomb, 1 Confined, 2 Higgs(Zq-deconf), 3 Higgs(deep Zq-conf)
def classify(r):
    higgs = r['cos'] > COS_ORD
    conf1 = (not np.isnan(r['sig1'])) and r['sig1'] > SIG_CONF
    if higgs:
        return 3 if conf1 else 2
    return 1 if conf1 else 0

for r in recs: r['ph'] = classify(r)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
except Exception as e:
    print("no matplotlib:", e); raise SystemExit
os.makedirs("u1f_campaign_analysis", exist_ok=True)
cmap = ListedColormap(["#3b6fb0", "#b03b3b", "#6fb03b", "#206020"])  # Coulomb, Confined, Higgs-deconf, Higgs-deepZq
labels = ["Coulomb (sigma1~0,|P1|>0)", "Confined (sigma1>0,cos low,rhoM hi)",
          "Higgs Zq-deconf (cos hi,sigma1~0)", "Higgs deep-Zq-conf (cos hi,sigma1>0,sigmaq~0)"]

qs = [2, 3, 4, 5, 6, 8]
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
for ax, q in zip(axes.flat, qs):
    pts = [r for r in recs if r['q'] == q and r['L'] == 8]
    if not pts:
        ax.set_title(f"q={q} (no data)"); continue
    bs = sorted(set(r['beta'] for r in pts)); ks = sorted(set(r['kappa'] for r in pts))
    grid = np.full((len(ks), len(bs)), np.nan)
    for r in pts: grid[ks.index(r['kappa']), bs.index(r['beta'])] = r['ph']
    ax.imshow(grid, origin="lower", aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5,
              extent=[min(bs), max(bs), min(ks), max(ks)])
    ax.set_xlabel("beta"); ax.set_ylabel("kappa"); ax.set_title(f"q={q}  ({len(pts)} pts)")
fig.legend(handles=[Patch(facecolor=cmap(i), label=labels[i]) for i in range(4)],
           loc="lower center", ncol=2, fontsize=9)
fig.suptitle("Frozen U(1)+charge-q Higgs phase map (L=8) -- DEFINITIVE discriminants", fontsize=13)
fig.tight_layout(rect=[0, 0.07, 1, 0.97])
p = "u1f_campaign_analysis/phasemap_definitive.png"; fig.savefig(p, dpi=110); print("wrote", p)

# triple-point-ish summary: for each q, list where Coulomb(0), Confined(1), Higgs(2/3) coexist as neighbors
print("\n## phase census per q (L=8): counts {Coulomb, Confined, Higgs-deconf, Higgs-deepZq}")
for q in qs:
    pts = [r for r in recs if r['q'] == q and r['L'] == 8]
    if not pts: continue
    c = [sum(1 for r in pts if r['ph'] == i) for i in range(4)]
    # Coulomb persistence: max kappa at which any Coulomb point exists (wedge => large)
    coul_k = [r['kappa'] for r in pts if r['ph'] == 0]
    kmax_coul = max(coul_k) if coul_k else None
    print(f"  q={q}: Coul={c[0]} Conf={c[1]} Higgs_dec={c[2]} Higgs_deepZq={c[3]}  | Coulomb persists to kappa<= {kmax_coul}")
