#!/usr/bin/env python3
"""Plot a fixed-kappa beta-scan (the freezing probe): chi_plaq(beta) and <plaq>(beta)
per rep, from the per-trajectory ts files, with blocking error bars. The freezing of
the residual discrete gauge theory should show as a chi_plaq (gauge specific-heat) bump.
Usage: python3 scripts/su2_freeze_plot.py <scan_dir> [kappa_tag]   (default su2scan_freeze, k12)"""
import sys, os, glob, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "su2scan_freeze"
KTAG = sys.argv[2] if len(sys.argv) > 2 else "12"
V = 16
reps = ["2T", "2O", "2I"]
betaf = {"2T": 2.24, "2O": 3.26, "2I": 5.82}   # literature freezing beta_f

def series(rep):
    out = []
    for f in glob.glob(os.path.join(D, "ts", f"{rep}_b*_k{KTAG}.dat")):
        m = re.search(rf"{rep}_b([\d.]+)_k{KTAG}\.dat$", os.path.basename(f))
        if not m: continue
        try:
            d = np.loadtxt(f, comments="#")
            if d.ndim != 2 or len(d) < 20: continue
            p = d[:, 1]
            chi = V * p.var()
            nb = 20; bs = len(p)//nb
            chie = np.std([V*np.var(p[i*bs:(i+1)*bs]) for i in range(nb)])/np.sqrt(nb) if bs > 5 else 0.0
            out.append((float(m.group(1)), p.mean(), chi, chie))
        except Exception: pass
    out.sort()
    return out

present = [(r, series(r)) for r in reps]
present = [(r, s) for r, s in present if s]
if not present:
    print(f"no ts data in {D}/ts for k{KTAG}"); sys.exit(1)

fig, ax = plt.subplots(2, len(present), figsize=(5*len(present), 7), squeeze=False)
for j, (r, s) in enumerate(present):
    b = np.array([x[0] for x in s]); pl = np.array([x[1] for x in s])
    chi = np.array([x[2] for x in s]); chie = np.array([x[3] for x in s])
    ax[0][j].errorbar(b, chi, yerr=chie, fmt="s-", color="crimson", capsize=3)
    ax[0][j].axvline(betaf.get(r, np.nan), ls="--", color="gray", label=f"expected freeze β_f={betaf.get(r)}")
    ax[0][j].set(title=f"{r}: χ_plaq(β) at κ={KTAG}", xlabel="β", ylabel="χ_plaq")
    ax[0][j].legend(fontsize=8); ax[0][j].grid(alpha=0.3)
    ax[1][j].plot(b, pl, "o-", color="navy")
    ax[1][j].axvline(betaf.get(r, np.nan), ls="--", color="gray")
    ax[1][j].set(title=f"{r}: ⟨plaq⟩(β) at κ={KTAG}", xlabel="β", ylabel="⟨plaq⟩")
    ax[1][j].grid(alpha=0.3)
fig.suptitle(f"Freezing probe: fixed κ={KTAG} β-scan (dashed = literature β_f); a freezing = χ_plaq bump")
fig.tight_layout(rect=[0,0,1,0.96])
out = os.path.join(D, f"freeze_k{KTAG}.png"); fig.savefig(out, dpi=120); print("wrote", out)
for r, s in present:
    print(f"\n{r}:  beta  <plaq>   chi_plaq +/- err")
    for b, pl, chi, chie in s: print(f"   {b:5.2f}  {pl:.3f}   {chi:.4f} +/- {chie:.4f}")
