#!/usr/bin/env python3
"""SU(2)->H 3-PHASE CLASSIFICATION map from the observables we already have, using cutoffs
(works for crossovers too -- a cutoff just picks a side). Each (beta,kappa) -> one of:
  HIGGS    : matter condensed, L_link > L_cut.
  COULOMB  : not Higgs, and beta > beta_c (above the gauge transition).
  CONFINED : not Higgs, and beta < beta_c.
beta_c = the gauge transition location = chi_plaq (gauge specific heat) PEAK ("chi location"),
taken as a single robust value = median of the per-kappa peaks over the non-Higgs rows.
Regions meet at a triple point by construction (no crossing lines).

MULTI-SOURCE: pass a comma-list of scan dirs; each rep is loaded from the FIRST dir that has
it well-populated (>= MINPTS points), else the next -- so e.g. adj/2T/2O come from a high-stat
dir and 2I falls back to the full-grid dir.
Usage: python3 scripts/su2_phase_classify.py [dir1,dir2,...] [L_cut=0.5]"""
import sys, os, csv, glob, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

DIRS = (sys.argv[1] if len(sys.argv) > 1 else "su2scan_v2_hot,su2scan_grid").split(",")
L_CUT = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
MINPTS = 40
V = 16
ORDER = [("adj", "adjoint -> U(1)"), ("2T", "2T (spin-3)"), ("2O", "2O (spin-4)"), ("2I", "2I (spin-6)")]

def chi_plaq_map(rep, D):
    out = {}
    for f in glob.glob(os.path.join(D, "ts", f"{rep}_*.dat")):
        try:
            hdr = open(f).readline()
            mb = re.search(r"beta=([\d.]+)", hdr); mk = re.search(r"kappa=([\d.]+)", hdr)
            d = np.loadtxt(f, comments="#")
            if mb and mk and d.ndim == 2 and len(d) >= 10:
                out[(round(float(mb.group(1)), 3), round(float(mk.group(1)), 3))] = V * np.var(d[:, 1])
        except Exception: pass
    return out

def load(rep, D):
    p = os.path.join(D, f"{rep}.csv")
    if not os.path.exists(p): return None
    cp = chi_plaq_map(rep, D); rows = []
    for r in csv.DictReader(open(p)):
        try:
            b, k = round(float(r["beta"]), 3), round(float(r["kappa"]), 3)
            rows.append({"beta": b, "kappa": k, "Llink": float(r["Llink"]),
                         "chi_plaq": cp.get((b, k), np.nan)})
        except Exception: pass
    return rows or None

def pick(rep):
    best = None
    for D in DIRS:
        rows = load(rep, D)
        if rows and len(rows) >= MINPTS: return D, rows           # first well-populated
        if rows and (best is None or len(rows) > len(best[1])): best = (D, rows)
    return best                                                    # else most-populated, or None

def classify(rows, L_cut):
    bs = sorted({r["beta"] for r in rows}); ks = sorted({r["kappa"] for r in rows})
    Ll = {(r["beta"], r["kappa"]): r["Llink"] for r in rows}
    cp = {(r["beta"], r["kappa"]): r["chi_plaq"] for r in rows}
    betac = {}
    for k in ks:
        prof = [(b, cp.get((b, k), np.nan)) for b in bs]; prof = [(b, v) for b, v in prof if not np.isnan(v)]
        if len(prof) >= 3 and max(v for _, v in prof) > 0:
            betac[k] = prof[int(np.argmax([v for _, v in prof]))][0]
    nonh = [k for k in ks if np.mean([Ll.get((b, k), 0.0) for b in bs]) < L_cut]
    src = [betac[k] for k in nonh if k in betac] or list(betac.values()) or [float(np.median(bs))]
    bc = float(np.median(src))
    P = np.full((len(ks), len(bs)), np.nan)
    for i, k in enumerate(ks):
        for j, b in enumerate(bs):
            L = Ll.get((b, k), np.nan)
            if np.isnan(L): continue
            P[i, j] = 2 if L > L_cut else (1 if b > bc else 0)
    return np.array(bs), np.array(ks), P, bc

present = []
for lab, pretty in ORDER:
    got = pick(lab)
    if got: present.append((lab, pretty, got[0], got[1]))
cmap = ListedColormap(["#3b6ea5", "#5aae61", "#c0392b"])
fig, axes = plt.subplots(1, len(present), figsize=(4.5 * len(present), 4.5), squeeze=False)
for j, (lab, pretty, src, rows) in enumerate(present):
    bs, ks, P, bc = classify(rows, L_CUT)
    ax = axes[0][j]
    ax.imshow(P, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=2,
              extent=[bs.min() - .25, bs.max() + .25, ks.min() - .5, ks.max() + .5])
    ax.axvline(bc, color="k", ls=":", lw=1, alpha=0.5)
    ax.set(title=f"{pretty}  (β_c≈{bc:.1f})\n[{src}]", xlabel="β", ylabel="κ")
    ax.text(0.22 * bs.max(), 0.10 * ks.max(), "confined", color="white", fontsize=9, ha="center")
    ax.text(0.82 * bs.max(), 0.10 * ks.max(), "Coulomb", color="white", fontsize=9, ha="center")
    ax.text(0.55 * bs.max(), 0.82 * ks.max(), "HIGGS", color="white", fontsize=12, ha="center", weight="bold")
fig.legend(handles=[Patch(color="#3b6ea5", label="confined"), Patch(color="#5aae61", label="Coulomb"),
                    Patch(color="#c0392b", label="Higgs")], loc="upper right", ncol=3, fontsize=9)
fig.suptitle(f"SU(2)→H 3-phase classification (L=2): Higgs if L_link>{L_CUT}; confined/Coulomb split at χ_plaq-peak β_c")
fig.tight_layout(rect=[0, 0, 1, 0.92])
out = "su2_phase_classify_v2.png"; fig.savefig(out, dpi=120); print("wrote", out)
for lab, pretty, src, rows in present:
    bs, ks, P, bc = classify(rows, L_CUT)
    frac = {n: int((P == v).sum()) for v, n in [(0, "conf"), (1, "Coul"), (2, "Higgs")]}
    print(f"  {lab:5s} <-{src:16s} npts={len(rows):3d} β_c≈{bc:.2f}  {frac}")
