#!/usr/bin/env python3
"""Plot the SU(2)->H (beta,kappa) wide scan: per-rep heatmaps of
  L_link (scalar order param) / chi_link (Higgs transition) / chi_plaq (gauge crossover) / plaq.
chi_plaq = V*Var(plaquette) is computed from the per-trajectory ts files (the gauge-sector
specific heat); the other three come from the per-rep CSV. White dashed = chi_link ridge.
Usage: python3 scripts/su2_phase_plot.py [scan_dir]   (default su2scan_w)"""
import sys, os, csv, glob, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

DIRN = sys.argv[1] if len(sys.argv) > 1 else "su2scan_w"
V = 16  # 2^4 sites (matches the driver's susceptibility normalization)
ORDER = [("adj", "adjoint -> U(1) [2^4]"), ("Q8soft", "spin-2 -> Q8 soft [2^4]"),
         ("2T", "2T (spin-3) [2^4]"), ("2O", "2O (spin-4) [2^4]"), ("2I", "2I (spin-6) [2^4]")]

def chi_plaq_map(label):
    """chi_plaq = V*Var(plaq) per (beta,kappa), read from ts headers + the plaq column."""
    out = {}
    for f in glob.glob(os.path.join(DIRN, "ts", f"{label}_*.dat")):
        try:
            with open(f) as fh: hdr = fh.readline()
            mb = re.search(r"beta=([\d.]+)", hdr); mk = re.search(r"kappa=([\d.]+)", hdr)
            if not (mb and mk): continue
            d = np.loadtxt(f, comments="#")
            if d.ndim != 2 or len(d) < 10: continue
            out[(round(float(mb.group(1)), 3), round(float(mk.group(1)), 3))] = V * np.var(d[:, 1])
        except Exception: pass
    return out

def load(label):
    path = os.path.join(DIRN, f"{label}.csv")
    if not os.path.exists(path): return None
    cp = chi_plaq_map(label)
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                row = {k: float(r[k]) for k in ("beta", "kappa", "plaq", "Llink", "chi_link")}
                row["chi_plaq"] = cp.get((round(row["beta"], 3), round(row["kappa"], 3)), np.nan)
                rows.append(row)
            except (ValueError, KeyError): pass
    return rows or None

def grid(rows, key):
    bs = sorted({r["beta"] for r in rows}); ks = sorted({r["kappa"] for r in rows})
    M = np.full((len(ks), len(bs)), np.nan)
    for r in rows: M[ks.index(r["kappa"]), bs.index(r["beta"])] = r[key]
    return np.array(bs), np.array(ks), M

def higgs_boundary(rows):
    """kappa_c(beta): kappa of the steepest rise of <L_link> (matter sector)."""
    out = []
    for b in sorted({r["beta"] for r in rows}):
        kv = sorted((r["kappa"], r["Llink"]) for r in rows if r["beta"] == b)
        if len(kv) < 2: continue
        ks = [k for k, _ in kv]; vs = [v for _, v in kv]
        sl = [(vs[i+1]-vs[i])/(ks[i+1]-ks[i]) for i in range(len(ks)-1)]
        i = int(np.argmax(sl)); out.append((b, 0.5*(ks[i]+ks[i+1])))
    return np.array(out) if out else np.empty((0, 2))

def coulomb_boundary(rows):
    """beta_c(kappa): beta of the chi_plaq peak = gauge specific heat (gauge sector)."""
    out = []
    for k in sorted({r["kappa"] for r in rows}):
        bv = sorted((r["beta"], r["chi_plaq"]) for r in rows if r["kappa"] == k and not np.isnan(r["chi_plaq"]))
        if len(bv) < 3: continue
        bs = [b for b, _ in bv]; vs = np.array([v for _, v in bv])
        if vs.max() <= 0: continue
        out.append((k, bs[int(np.argmax(vs))]))
    return np.array(out) if out else np.empty((0, 2))

reps = [(l, p, load(l)) for l, p in ORDER]; reps = [(l, p, r) for l, p, r in reps if r]
if not reps: print(f"no CSVs in {DIRN}/"); sys.exit(1)
cols = [("Llink", "L_link (order param)", "viridis"),
        ("chi_link", "chi_link (Higgs transition)", "magma"),
        ("chi_plaq", "chi_plaq (gauge crossover)", "inferno"),
        ("plaq", "plaquette", "cividis")]
nr, nc = len(reps), len(cols)
fig, axes = plt.subplots(nr, nc, figsize=(4.1*nc, 3.0*nr), squeeze=False)
h_handle = c_handle = None
for i, (lab, pretty, rows) in enumerate(reps):
    hb = higgs_boundary(rows)     # (beta, kappa_c) -- matter/Higgs boundary
    cb = coulomb_boundary(rows)   # (kappa, beta_c) -- gauge/confining<->Coulomb boundary
    for j, (key, title, cmap) in enumerate(cols):
        ax = axes[i][j]; bs, ks, M = grid(rows, key)
        im = ax.imshow(M, origin="lower", aspect="auto", cmap=cmap,
                       extent=[bs.min()-.12, bs.max()+.12, ks.min()-.05, ks.max()+.05])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if len(hb): hl, = ax.plot(hb[:, 0], hb[:, 1], "o-", color="red", ms=3.5, lw=1.6, alpha=0.95); h_handle = hl
        if len(cb): cl, = ax.plot(cb[:, 1], cb[:, 0], "s--", color="cyan", ms=3.5, lw=1.6, alpha=0.95); c_handle = cl
        ax.set_title(f"{pretty}\n{title}", fontsize=8)
        if j == 0: ax.set_ylabel("kappa")
        if i == nr-1: ax.set_xlabel("beta")
if h_handle and c_handle:
    fig.legend([h_handle, c_handle],
               ["Higgs κ_c(β)  [L_link rise]", "confining↔Coulomb β_c(κ)  [χ_plaq peak]"],
               loc="upper right", fontsize=9, ncol=2)
fig.suptitle("SU(2)->H wide 2^4 scan: L_link / chi_link / chi_plaq / plaq   (red = Higgs boundary, cyan = confining↔Coulomb boundary)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(DIRN, "su2_phase.png"); fig.savefig(out, dpi=110); print("wrote", out)
print("\nboundaries:")
for lab, pretty, rows in reps:
    hb = higgs_boundary(rows); cb = coulomb_boundary(rows)
    print(f"  {lab:7s} Higgs κ_c(β)= " + " ".join(f"{b:.0f}:{k:.1f}" for b, k in hb))
    print(f"  {' '*7} conf↔Coul β_c(κ)= " + " ".join(f"{k:.0f}:{b:.0f}" for k, b in cb))
