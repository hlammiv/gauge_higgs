#!/usr/bin/env python3
"""Reconstruct rho(A,B) from a RECTANGULAR 2D-LLR tile (u1_llr rect) and map the
(beta,kappa) phase diagram + triple point.

Builds the (i,j) grid directly from the unique A and E2=-B cell-centre values (the rect
raster is a regular grid), integrates ln rho via the same weighted-least-squares gradient
solve as the ridge reconstruct, then reweights to a (beta,kappa) grid:
    <O>(beta,kappa) = sum_cells O * exp(lnrho - beta*A - kappa*E2) / Z      (E2=-B)
using the continuous-tile cell_moments (handles wide cells). Susceptibilities
chi_A=Var(A), chi_B=Var(B) peak along the gauge and Higgs transition lines; their
crossing is the triple point.

Usage: u1_llr_rect_reconstruct.py <rect.out> --q 2 --beta-range 0.5 2.0 --kappa-range 0 0.5
       [--heatmap out.png] [--amax 80]
"""
import sys, os, argparse, importlib.util
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("rec", os.path.join(_here, "u1_llr_reconstruct.py"))
rec = importlib.util.module_from_spec(spec); spec.loader.exec_module(rec)


def parse_rect(fname, amax):
    cells = []; geo = {}
    for l in open(fname):
        if l.startswith("LLR2DRECT"):
            for tok in l.replace("|", " ").split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    try: geo[k] = float(v)
                    except ValueError: pass
            import re
            for key, pat in (("hw1", r"hw1=([-\d.eE+]+)"), ("hw2", r"hw2=([-\d.eE+]+)")):
                m = re.search(pat, l)
                if m: geo[key] = float(m.group(1))
        elif l.startswith("ANE2:"):
            t = l.split(); a1, a2 = float(t[3]), float(t[4])
            if np.isfinite(a1) and np.isfinite(a2) and abs(a1) < amax and abs(a2) < amax:
                cells.append([float(t[1]), float(t[2]), a1, a2])
    NP, NL, V, D = rec.norms_from_header([fname])  # falls back if header differs
    return np.array(cells), geo


def build_grid(cells):
    """(i,j) from unique sorted A (descending -> i) and E2 (ascending -> j)."""
    A = cells[:, 0]; E2 = cells[:, 1]
    Au = np.unique(np.round(A, 3))[::-1]      # high A -> i=0
    Eu = np.unique(np.round(E2, 3))           # low E2 (high B) -> j=0
    ai = {a: i for i, a in enumerate(Au)}
    ej = {e: j for j, e in enumerate(Eu)}
    idx = [(ai[round(a, 3)], ej[round(e, 3)]) for a, e in zip(A, E2)]
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--q", type=int, default=2)
    ap.add_argument("--amax", type=float, default=80.0)
    ap.add_argument("--beta-range", type=float, nargs=2, default=[0.5, 2.0])
    ap.add_argument("--kappa-range", type=float, nargs=2, default=[0.0, 0.5])
    ap.add_argument("--ngrid", type=int, default=60)
    ap.add_argument("--e2-weight", type=float, default=0.3,
                    help="downweight a2 (E2-direction) edges in the lnrho integration; <1 keeps "
                         "the gauge (A) direction governed by the clean a1 so beta_c lands right")
    ap.add_argument("--heatmap", default=None)
    args = ap.parse_args()

    cells, geo = parse_rect(args.file, args.amax)
    if len(cells) < 6:
        sys.exit(f"FAIL: only {len(cells)} cells")
    idx = build_grid(cells)
    lnrho = rec.integrate(cells, idx, e2_weight=args.e2_weight)
    A = cells[:, 0]; E2 = cells[:, 1]
    hw1 = geo.get("hw1", 1.0); hw2 = geo.get("hw2", 1.0)
    print(f"# rect: {len(cells)} cells, A in [{A.min():.0f},{A.max():.0f}], B in [{-E2.max():.0f},{-E2.min():.0f}], "
          f"a1 in [{cells[:,2].min():.2f},{cells[:,2].max():.2f}], a2 in [{cells[:,3].min():.3f},{cells[:,3].max():.3f}]")

    blo, bhi = args.beta_range; klo, khi = args.kappa_range
    betas = np.linspace(blo, bhi, args.ngrid); kappas = np.linspace(klo, khi, args.ngrid)
    meanA = np.zeros((args.ngrid, args.ngrid)); meanB = np.zeros_like(meanA)
    chiA = np.zeros_like(meanA); chiB = np.zeros_like(meanA)
    for ik, k in enumerate(kappas):
        for ib, b in enumerate(betas):
            # physical exponent lnrho - beta*A - kappa*E2 (E2=-B): b1=beta, b2=+kappa
            S1, V1, S2, V2 = rec.cell_moments(cells, lnrho, hw1, hw2, b, k)
            meanA[ik, ib] = S1; meanB[ik, ib] = -S2; chiA[ik, ib] = V1; chiB[ik, ib] = V2

    # --- GAUGE line (confined<->Coulomb): extract from the LOW-B (symmetric) cells via a1 ONLY.
    # At low kappa B~0 on both sides (Clausius-Clapeyron: DeltaB~0 -> VERTICAL line), and a1 is the
    # clean gauge slope (the full 2D marginal mixes noisy-a2 B-rows -> mislocates it). beta_c =
    # steepest <A>(beta) of the 1D rho_gauge(A)=integral a1 dA along the lowest-B row.
    Bcells = -E2
    blo_mask = Bcells < (Bcells.min() + 0.18 * (Bcells.max() - Bcells.min()))   # symmetric/low-B slice
    As = A[blo_mask]; a1s = cells[blo_mask, 2]
    o = np.argsort(As); As = As[o]; a1s = a1s[o]
    # dedup repeated A (average a1)
    Au = np.unique(As); a1u = np.array([a1s[As == a].mean() for a in Au])
    lr = np.zeros_like(Au); lr[1:] = np.cumsum(0.5 * (a1u[1:] + a1u[:-1]) * np.diff(Au))
    mA = np.array([(Au * (lambda w: (e := np.exp(w - w.max())) / e.sum())(lr - b * Au)).sum() for b in betas])
    beta_c_gauge = float(betas[np.argmin(np.gradient(mA, betas))])
    print(f"\n# GAUGE line (confined<->Coulomb) from low-B a1-slice: beta_c = {beta_c_gauge:.3f}  "
          f"(VERTICAL; ClCl DeltaB~0).  [pure compact U(1) bulk ~1.01]")

    # --- HIGGS line kappa_c(beta) from steepest d<B>/dkappa (interior) ---
    edge = max(2, args.ngrid // 20)
    dBdk = np.abs(np.gradient(meanB, kappas, axis=0)); dBdk[:edge, :] = 0; dBdk[-edge:, :] = 0
    kc_higgs = kappas[np.argmax(dBdk, axis=0)]
    print("# HIGGS line kappa_c(beta) [steepest d<B>/dkappa]:")
    for ib in range(0, args.ngrid, max(1, args.ngrid // 10)):
        print(f"   beta={betas[ib]:.3f}  kappa_c={kc_higgs[ib]:.3f}")
    # triple point: where the vertical gauge line (beta=beta_c_gauge) meets the Higgs line
    ibg = int(np.argmin(np.abs(betas - beta_c_gauge)))
    kappa_t = float(kc_higgs[ibg])
    best = (beta_c_gauge, kappa_t)
    print(f"\n# TRIPLE-POINT estimate (vertical gauge line meets Higgs line): "
          f"beta_t~{beta_c_gauge:.3f}, kappa_t~{kappa_t:.3f}")
    bc_gauge = np.full_like(kappas, beta_c_gauge)      # vertical line for plotting

    if args.heatmap:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        NP = geo.get("NPLAQ", 1.0)
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
        kt = best[1]
        for a, Z, ttl in ((ax[0], meanA, "<A> (gauge)"), (ax[1], meanB, "<B> (matter)"),
                          (ax[2], np.log1p(chiB), "log chi_B (Higgs line)")):
            im = a.pcolormesh(betas, kappas, Z, shading="auto", cmap="viridis")
            # gauge line: vertical at beta_c, only BELOW the triple point (above it is Higgs)
            a.plot([beta_c_gauge, beta_c_gauge], [kappas.min(), kt], "w--", lw=1.5, alpha=0.85)
            # Higgs line: kappa_c(beta), only for beta >= triple point (left of it is the gauge line)
            mH = betas >= beta_c_gauge
            a.plot(betas[mH], kc_higgs[mH], "w:", lw=1.5, alpha=0.85)
            a.set_xlabel(r"$\beta$"); a.set_ylabel(r"$\kappa$"); a.set_title(ttl); fig.colorbar(im, ax=a)
        if best:
            for a in ax: a.plot([best[0]], [best[1]], "r*", ms=16)
        fig.suptitle(f"U(1)+charge-{args.q} phase diagram from rectangular 2D-LLR rho(A,B)")
        fig.tight_layout(); fig.savefig(args.heatmap, dpi=130, bbox_inches="tight")
        print("wrote", args.heatmap)


if __name__ == "__main__":
    main()
