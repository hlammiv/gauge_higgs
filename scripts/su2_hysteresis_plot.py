#!/usr/bin/env python3
"""Hysteresis-loop plot: overlay the HOT-start and COLD-start <plaq>(beta) branches at
fixed kappa for each rep. At L=2 high-kappa the HMC is basin-trapped, so a single-start
<plaq>(beta) is NOT an equilibrium curve -- but the SEPARATION between the hot and cold
branches IS an order-parameter-independent signal: where they split = metastable/coexistence
window (brackets a first-order freezing); where they merge = a single phase. The per-branch
error bar is the WITHIN-BASIN blocking error (valid for that branch's conditional mean only).
Usage: python3 scripts/su2_hysteresis_plot.py [hot_dir] [cold_dir] [kappa_tag]
       (defaults: su2scan_freeze  su2scan_hyst_cold  12)"""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HOT  = sys.argv[1] if len(sys.argv) > 1 else "su2scan_freeze"
COLD = sys.argv[2] if len(sys.argv) > 2 else "su2scan_hyst_cold"
KTAG = sys.argv[3] if len(sys.argv) > 3 else "12"
reps  = ["2T", "2O", "2I"]
betaf = {"2T": 2.24, "2O": 3.26, "2I": 5.82}   # literature freezing beta_f

def load(d, rep):
    f = os.path.join(d, f"{rep}.csv")
    if not os.path.isfile(f): return None
    try:
        a = np.genfromtxt(f, delimiter=",", names=True)
        a = np.atleast_1d(a)
        o = np.argsort(a["beta"])
        return a["beta"][o], a["plaq"][o], a["plaq_err"][o]
    except Exception:
        return None

present = []
for r in reps:
    h, c = load(HOT, r), load(COLD, r)
    if h is not None or c is not None:
        present.append((r, h, c))
if not present:
    print(f"no CSV data in {HOT} / {COLD}"); sys.exit(1)

fig, ax = plt.subplots(1, len(present), figsize=(5.2*len(present), 4.4), squeeze=False)
for j, (r, h, c) in enumerate(present):
    a = ax[0][j]
    if h is not None:
        a.errorbar(h[0], h[1], yerr=h[2], fmt="o-", color="crimson", capsize=3, label="hot start (disordered)")
    if c is not None:
        a.errorbar(c[0], c[1], yerr=c[2], fmt="s-", color="navy", capsize=3, label="cold start (ordered)")
    a.axvline(betaf.get(r, np.nan), ls="--", color="gray", label=f"lit. β_f={betaf.get(r)}")
    a.set(title=f"{r}: hysteresis loop at κ={KTAG}", xlabel="β", ylabel="⟨plaq⟩", ylim=(0, 1))
    a.legend(fontsize=8); a.grid(alpha=0.3)
fig.suptitle("SU(2)→H freezing: hot/cold ⟨plaq⟩ branches (gap = first-order coexistence window)")
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = f"su2_hysteresis_k{KTAG}.png"; fig.savefig(out, dpi=120); print("wrote", out)

for r, h, c in present:
    print(f"\n{r}:  beta   plaq_hot   plaq_cold   gap")
    bs = sorted(set((h[0] if h is not None else []).tolist()) | set((c[0] if c is not None else []).tolist()))
    hd = dict(zip(h[0], h[1])) if h is not None else {}
    cd = dict(zip(c[0], c[1])) if c is not None else {}
    for b in bs:
        ph, pc = hd.get(b), cd.get(b)
        gap = f"{abs(pc-ph):.3f}" if (ph is not None and pc is not None) else "  -"
        print(f"   {b:5.2f}  {ph if ph is None else f'{ph:.4f}':>8}   {pc if pc is None else f'{pc:.4f}':>8}   {gap:>6}")
