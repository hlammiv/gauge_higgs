#!/usr/bin/env python3
"""High-kappa hot/cold hysteresis: at fixed high kappa, scan beta from a HOT (disordered)
and a COLD (ordered) start. Where the two <plaq>(beta) branches SPLIT is the metastable
coexistence window that brackets the confining<->Coulomb(freezing) gauge transition beta_c
-- the locator the single-start chi_plaq peak couldn't pin at high kappa.
Rows = reps, cols = high kappa. Usage: python3 scripts/su2_hyst_highk.py [hot_dir] [cold_dir]"""
import sys, os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HOT  = sys.argv[1] if len(sys.argv) > 1 else "su2scan_v2_hot"
COLD = sys.argv[2] if len(sys.argv) > 2 else "su2scan_v2_cold"
reps  = ["adj", "2T", "2O", "2I"]
betaf = {"2T": 2.24, "2O": 3.26, "2I": 5.82}
KS = [8.0, 12.0, 20.0]

def load(d, rep):
    p = os.path.join(d, f"{rep}.csv"); out = {}
    if not os.path.exists(p): return out
    for r in csv.DictReader(open(p)):
        try: out[(round(float(r["beta"]),3), round(float(r["kappa"]),3))] = (float(r["plaq"]), float(r["plaq_err"]))
        except Exception: pass
    return out

present = [r for r in reps if load(HOT, r) and load(COLD, r)]
fig, axes = plt.subplots(len(present), len(KS), figsize=(4.2*len(KS), 3.1*len(present)), squeeze=False)
for i, rep in enumerate(present):
    H = load(HOT, rep); C = load(COLD, rep)
    for j, k in enumerate(KS):
        ax = axes[i][j]
        hb = sorted((b, *H[(b,kk)]) for (b,kk) in H if kk == k)
        cb = sorted((b, *C[(b,kk)]) for (b,kk) in C if kk == k)
        if hb: b,p,e = zip(*hb); ax.errorbar(b,p,yerr=e,fmt="o-",color="crimson",ms=3,lw=1.4,capsize=2,label="hot start")
        if cb: b,p,e = zip(*cb); ax.errorbar(b,p,yerr=e,fmt="s--",color="navy",ms=3,lw=1.4,capsize=2,label="cold start")
        if rep in betaf: ax.axvline(betaf[rep], ls=":", color="gray", alpha=.8, label=f"lit. β_f={betaf[rep]}")
        ax.set(title=f"{rep}  κ={k:.0f}", ylim=(0,1)); ax.grid(alpha=0.3)
        if i == len(present)-1: ax.set_xlabel("β")
        if j == 0: ax.set_ylabel("⟨plaq⟩")
        if i == 0 and j == 0: ax.legend(fontsize=7)
fig.suptitle("High-κ hot/cold hysteresis: hot & cold ⟨plaq⟩ split = confining↔Coulomb (freezing) coexistence bracket")
fig.tight_layout(rect=[0,0,1,0.96])
out = os.path.join(HOT, "su2_hyst_highk.png"); fig.savefig(out, dpi=120); print("wrote", out)
for rep in present:
    H = load(HOT, rep); C = load(COLD, rep)
    print(f"\n{rep}: confining↔Coulomb bracket (β where hot/cold ⟨plaq⟩ gap > 0.05):")
    for k in KS:
        bs = sorted({b for b,kk in H if kk == k} & {b for b,kk in C if kk == k})
        sp = [b for b in bs if abs(H[(b,k)][0]-C[(b,k)][0]) > 0.05]
        print(f"   κ={k:.0f}: split at β= {sp if sp else '(none - branches agree, single phase)'}")
