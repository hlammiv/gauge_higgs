#!/usr/bin/env python3
"""SU(2) kappa=0 crossover from existing ts files (NO rerun). At kappa=0 the Higgs
decouples, so every rep's kappa=0 gauge ensemble is the SAME pure-SU(2) Wilson theory
(independent seeds) -> pool them for stats. Plot <plaq>(beta) [smooth] and the plaquette
specific heat chi_plaq = V*Var(plaq) [broad crossover bump, NOT a divergent transition].
Usage: python3 scripts/su2_crossover.py [scan_dir]"""
import sys, glob, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else "su2scan_w"
V = 16  # 2^4 sites (matches the driver's susceptibility normalization)
betas = ["1.00","1.25","1.50","1.75","2.00","2.25","2.50","2.75","3.00"]
reps = ["adj","Q8soft","2T","2O","2I"]

B, PL, CHI, CHIe = [], [], [], []
for b in betas:
    pooled = []
    for r in reps:
        f = os.path.join(D, "ts", f"{r}_b{b}_k0.0.dat")
        if os.path.exists(f):
            try:
                d = np.loadtxt(f, comments="#")
                if d.ndim == 2 and len(d) > 10: pooled.append(d[:, 1])  # plaq column
            except Exception: pass
    if not pooled: continue
    x = np.concatenate(pooled)
    B.append(float(b)); PL.append(x.mean()); CHI.append(V * x.var())
    # blocking error on chi_plaq: split into 20 blocks, std of per-block V*Var
    nb = 20; bs = len(x)//nb
    if bs > 5:
        chis = [V*np.var(x[i*bs:(i+1)*bs]) for i in range(nb)]
        CHIe.append(np.std(chis)/np.sqrt(nb))
    else: CHIe.append(0.0)

B = np.array(B)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(B, PL, "o-", color="navy")
ax[0].set(xlabel=r"$\beta$", ylabel=r"$\langle$plaq$\rangle$",
          title=r"SU(2) $\kappa{=}0$: $\langle$plaq$\rangle$ rises smoothly (no jump)")
ax[0].grid(alpha=0.3)
ax[1].errorbar(B, CHI, yerr=CHIe, fmt="s-", color="crimson", capsize=3)
ax[1].set(xlabel=r"$\beta$", ylabel=r"$\chi_{\rm plaq}=V\,{\rm Var}(\rm plaq)$",
          title=r"SU(2) $\kappa{=}0$: specific heat — broad crossover bump, no divergence")
ax[1].grid(alpha=0.3)
fig.suptitle("Pure SU(2) gauge crossover (2$^4$, kappa=0, all reps pooled): a CROSSOVER, not a confinement transition",
             fontsize=10)
fig.tight_layout(rect=[0,0,1,0.95])
out = os.path.join(D, "su2_crossover_k0.png")
fig.savefig(out, dpi=120); print("wrote", out)
print("beta   <plaq>   chi_plaq")
for i in range(len(B)): print(f"  {B[i]:.2f}  {PL[i]:.4f}  {CHI[i]:.4f}")
