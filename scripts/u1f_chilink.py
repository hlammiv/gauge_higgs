#!/usr/bin/env python3
"""Locate the Higgs transition kappa_c(q) from the matter susceptibility chi_link = Var(B)/vol peak
(the right order-parameter observable; <cos>_link is the q-blind hopping ENERGY). Reads the high-kappa
extension (u1f_ext + u1f_ext_lenore)."""
import glob, re, os
import numpy as np
HDR = re.compile(r"L=(\d+)\s+beta=([\d.]+)\s+kappa=([\d.]+)\s+q=(\d+)")
CHI = re.compile(r"chi_plaq=([\-\d.]+)\s+chi_link=([\-\d.]+)")
COS = re.compile(r"<cos>_link=([\-\d.]+)")

recs = []
for d in ("u1f_ext", "u1f_ext_lenore"):
    for f in glob.glob(d + "/job_point_*.out"):
        t = open(f).read(); h, c = HDR.search(t), CHI.search(t)
        if h and c:
            cs = COS.search(t)
            recs.append(dict(L=int(h[1]), beta=float(h[2]), kappa=float(h[3]), q=int(h[4]),
                             chiP=float(c[1]), chiL=float(c[2]), cos=float(cs[1]) if cs else np.nan))
print(f"# extension points with chi: {len(recs)}")

def peak_kappa(rows):
    rows = sorted(rows, key=lambda r: r['kappa'])
    if len(rows) < 3: return None, None
    ks = np.array([r['kappa'] for r in rows]); ch = np.array([r['chiL'] for r in rows])
    i = int(np.argmax(ch))
    return ks[i], ch[i]

print("\n## kappa_c(q) from chi_link PEAK  (Higgs transition; does it grow with q?)")
for b in (0.60, 0.90, 1.10):
    line = f"  beta={b}:  "
    for q in (2, 4, 6, 8):
        rows = [r for r in recs if r['q'] == q and r['L'] == 8 and abs(r['beta'] - b) < 1e-6]
        kc, ch = peak_kappa(rows)
        line += f"q{q}:kc={kc if kc is not None else '--':<5} " if kc is None else f"q{q}:kc={kc:<4g}(chi={ch:.1f})  "
    print(line)

print("\n## full chi_link(kappa) at beta=0.90, L=8 (peak = transition):")
for q in (2, 4, 6, 8):
    rows = sorted([r for r in recs if r['q'] == q and r['L'] == 8 and abs(r['beta'] - 0.90) < 1e-6], key=lambda r: r['kappa'])
    print(f"  q={q}: " + " ".join(f"{r['kappa']:g}:{r['chiL']:.1f}" for r in rows))

# plot chi_link(kappa) per q at beta=0.9
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    os.makedirs("u1f_campaign_analysis", exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for q in (2, 4, 6, 8):
        rows = sorted([r for r in recs if r['q'] == q and r['L'] == 8 and abs(r['beta'] - 0.90) < 1e-6], key=lambda r: r['kappa'])
        if not rows: continue
        ks = [r['kappa'] for r in rows]
        ax[0].plot(ks, [r['chiL'] for r in rows], 'o-', label=f"q={q}")
        ax[1].plot(ks, [r['cos'] for r in rows], 'o-', label=f"q={q}")
    ax[0].set_xlabel("kappa"); ax[0].set_ylabel("chi_link = Var(B)/vol"); ax[0].set_title("matter SUSCEPTIBILITY (peak=Higgs transition)"); ax[0].set_xscale('log'); ax[0].legend()
    ax[1].set_xlabel("kappa"); ax[1].set_ylabel("<cos>_link"); ax[1].set_title("hopping ENERGY (q-blind)"); ax[1].set_xscale('log'); ax[1].legend()
    fig.suptitle("frozen U(1)+charge-q Higgs, beta=0.9, L=8: susceptibility vs energy")
    fig.tight_layout(); p = "u1f_campaign_analysis/chilink_vs_kappa_b0.9.png"; fig.savefig(p, dpi=120); print("\nwrote", p)
except Exception as e:
    print("plot skipped:", e)
