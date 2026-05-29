#!/usr/bin/env python3
"""Analyze the SU(2)-fundamental Higgs finite-size scan (validation_su2higgs.dat).

Reports, per volume: the L_link curve, the susceptibility chi_link peak location
(the transition estimate), and saves a plot. The chi_link peak growing/sharpening
with volume is the finite-size signature of the confinement<->Higgs transition.
"""
import sys
from collections import defaultdict

fn = sys.argv[1] if len(sys.argv) > 1 else "validation_su2higgs.dat"
rows = defaultdict(list)  # L -> list of (kappa, L_link, err, chi_link, chi_phi, binder, acc)
with open(fn) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        try:
            L = int(p[0]); k = float(p[1]); ll = float(p[2]); er = float(p[3])
            cl = float(p[4]); cp = float(p[5]); bn = float(p[6]); ac = float(p[7])
        except (ValueError, IndexError):
            continue
        rows[L].append((k, ll, er, cl, cp, bn, ac))

print(f"{'L':>3} {'kappa':>7} {'L_link':>9} {'err':>8} {'chi_link':>9} {'chi_phi':>9} {'Binder':>7} {'acc':>5}")
for L in sorted(rows):
    for (k, ll, er, cl, cp, bn, ac) in sorted(rows[L]):
        mark = ""
        print(f"{L:>3} {k:>7.3f} {ll:>9.4f} {er:>8.4f} {cl:>9.4f} {cp:>9.4f} {bn:>7.3f} {ac:>5.2f}{mark}")
    print()

print("=== transition estimate (chi_link peak) per volume ===")
for L in sorted(rows):
    data = sorted(rows[L])
    kpk, clpk = max(((k, cl) for (k, ll, er, cl, cp, bn, ac) in data), key=lambda t: t[1])
    print(f"  L={L}: chi_link peaks at kappa ~ {kpk:.3f}  (chi_link = {clpk:.3f})")

# Optional plot
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for L in sorted(rows):
        d = sorted(rows[L])
        ks = [r[0] for r in d]
        ax[0].errorbar(ks, [r[1] for r in d], yerr=[r[2] for r in d], marker="o", label=f"$L={L}$")
        ax[1].plot(ks, [r[3] for r in d], marker="s", label=f"$L={L}$")
    ax[0].set_xlabel(r"$\kappa$"); ax[0].set_ylabel(r"$L_{\rm link}$"); ax[0].set_title("gauge-invariant link energy"); ax[0].legend()
    ax[1].set_xlabel(r"$\kappa$"); ax[1].set_ylabel(r"$\chi_{\rm link}=V\,{\rm Var}(L_{\rm link})$"); ax[1].set_title("susceptibility (peak = transition)"); ax[1].legend()
    fig.suptitle(r"SU(2)-fundamental Higgs, $\beta=2.3$, $\lambda=0.5$ — finite-size scan")
    fig.tight_layout()
    fig.savefig("validation_su2higgs.png", dpi=120)
    print("\nplot -> validation_su2higgs.png")
except Exception as e:
    print(f"\n(plot skipped: {e})")
