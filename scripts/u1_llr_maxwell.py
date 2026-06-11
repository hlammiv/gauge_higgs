#!/usr/bin/env python3
"""Maxwell equal-weight + latent-heat + free-energy-barrier analysis of an LLR rho(A,B) tile.

NUMERICALLY OBSERVES a transition instead of drawing a line: builds the reweighted
probability of an order parameter, P(O) propto sum_cells exp(lnrho - beta*A + kappa*B)
binned in O, finds its TWO peaks, and locates the coupling where they have EQUAL WEIGHT
(Maxwell). Reports:
  - transition coupling (equal-weight),
  - latent heat = separation of the two peaks in O,
  - barrier  Delta F = peak - intervening valley of -ln P  (the quantity that GROWS with
    volume for a genuine first-order transition; flat/shrinking => crossover).
Run per-L and compare Delta F(L) for the finite-size-scaling order test.

  higgs mode: order param = B, fix beta, scan kappa.
  gauge mode: order param = A, fix kappa, scan beta.

Usage: u1_llr_maxwell.py <combined.out> --mode higgs --fix-beta 1.0 --scan 0.2 0.6
                          [--plot out.png] [--amax 80]
"""
import sys, os, argparse, importlib.util
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
def _load(name, fn):
    s = importlib.util.spec_from_file_location(name, os.path.join(_here, fn))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
rec = _load("rec", "u1_llr_reconstruct.py")
rr  = _load("rr", "u1_llr_rect_reconstruct.py")


def logsumexp_bins(vals, keys, ukeys):
    out = np.full(len(ukeys), -np.inf)
    for i, u in enumerate(ukeys):
        v = vals[keys == u]
        if len(v):
            m = v.max(); out[i] = m + np.log(np.exp(v - m).sum())
    return out


def two_peaks(P):
    """Two highest INTERIOR local maxima separated by a valley. ENDPOINTS ARE NOT PEAKS --
    a maximum at the domain edge means the tile is too small (the real peak is outside), so it
    is a boundary artifact, not a measurement. Returns ((i1,i2,iv) or None, edge_pileup_flag)."""
    n = len(P); fin = np.isfinite(P)
    edge = (fin[0] and fin[1] and P[0] >= P[1]) or (fin[n-1] and fin[n-2] and P[n-1] >= P[n-2])
    pk = [i for i in range(1, n - 1) if fin[i] and P[i] >= P[i-1] and P[i] >= P[i+1]]
    if len(pk) < 2:
        return None, edge
    pk.sort(key=lambda i: -P[i]); i1, i2 = sorted(pk[:2])
    iv = i1 + int(np.argmin(P[i1:i2+1]))
    if iv == i1 or iv == i2:
        return None, edge
    return (i1, i2, iv), edge


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--mode", choices=["higgs", "gauge"], default="higgs")
    ap.add_argument("--fix-beta", type=float, default=1.0)
    ap.add_argument("--fix-kappa", type=float, default=0.0)
    ap.add_argument("--scan", type=float, nargs=2, default=[0.2, 0.6])
    ap.add_argument("--nscan", type=int, default=121)
    ap.add_argument("--e2-weight", type=float, default=1.0)
    ap.add_argument("--amax", type=float, default=80.0)
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    cells, geo = rr.parse_rect(args.file, args.amax)
    idx = rr.build_grid(cells)
    lnrho = rec.integrate(cells, idx, e2_weight=args.e2_weight)
    A = cells[:, 0]; B = -cells[:, 1]
    NP, NL, V, D = rec.norms_from_header([args.file])

    if args.mode == "higgs":
        O = B; Ok = np.round(B, 0); uO = np.unique(Ok); fixed = f"beta={args.fix_beta}"
        def expo(s):  # s = kappa
            return lnrho - args.fix_beta * A + s * B
    else:
        O = A; Ok = np.round(A, 0); uO = np.unique(Ok); fixed = f"kappa={args.fix_kappa}"
        def expo(s):  # s = beta
            return lnrho - s * A + args.fix_kappa * B

    scans = np.linspace(args.scan[0], args.scan[1], args.nscan)
    dscan = (args.scan[1] - args.scan[0]) / (args.nscan - 1)
    spacing = float(np.median(np.diff(uO)))            # order-param grid quantization
    sym = "kappa_c" if args.mode == "higgs" else "beta_c"
    opar = "B" if args.mode == "higgs" else "A"
    best = None; n_dp = 0; n_edge = 0
    for s in scans:
        P = logsumexp_bins(expo(s), Ok, uO)
        tp, edge = two_peaks(P)
        if edge:
            n_edge += 1
        if tp is None:
            continue
        n_dp += 1
        i1, i2, iv = tp
        dh = P[i1] - P[i2]; barrier = min(P[i1], P[i2]) - P[iv]
        if best is None or abs(dh) < abs(best[1]):
            best = (s, dh, barrier, uO[i1], uO[i2], P, edge)

    print(f"# {args.mode} transition, fixed {fixed}, scan [{args.scan[0]},{args.scan[1]}] ({args.nscan} pts, d={dscan:.3f}), "
          f"{len(cells)} cells, V={V}")
    print(f"# order-param grid spacing d{opar}={spacing:.0f} -> peak positions/latent heat resolved only to +-{spacing:.0f}")
    print(f"# interior double-peak in {n_dp}/{args.nscan} scan pts; boundary pile-up in {n_edge}/{args.nscan}")
    if n_edge:
        print(f"# WARNING: probability reaches a domain edge in {n_edge} scan pts -> tile too small there; not physical.")
    if best is None:
        print(f"# RESULT: no interior double-peak in range -> transition NOT numerically resolved at this L/domain.")
        return
    s_c, dh, barrier, olo, ohi, P, edge = best
    balanced = abs(dh) < 2.0                            # equal-weight requires |dh|->0 (ln units)
    ok = balanced and not edge
    print(f"# best-balanced: {sym}={s_c:.3f}+-{dscan:.3f}(scan)  |peak-height diff|={abs(dh):.2f}"
          + ("" if balanced else "  <-- NOT balanced (|dh|>2): equal-weight NOT located"))
    print(f"# peaks at {opar}={olo:.0f},{ohi:.0f} (+-{spacing:.0f}); latent Delta{opar}={ohi-olo:.0f}+-{spacing:.0f}"
          + ("  [HIGH PEAK ON BOUNDARY -> artifact, not physical]" if edge else ""))
    print(f"# barrier Delta F={barrier:.2f} (ln P) at THIS SINGLE volume V={V}. "
          f"ORDER IS NOT DETERMINED FROM ONE L: first-order requires Delta F(L) to GROW with volume (FSS).")
    print(f"# VERDICT: {'genuine interior double-peak resolved' if ok else 'NOT a clean resolved transition (see warnings above)'}.")

    if args.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.figure(figsize=(6.5, 4.5)); Pn = P - P.max()
        plt.plot(uO, Pn, "o-", ms=3)
        plt.axvline(olo, color="g", ls=":"); plt.axvline(ohi, color="r", ls=":")
        plt.xlabel(f"{opar} ({'matter' if args.mode=='higgs' else 'gauge'})"); plt.ylabel("ln P (reweighted)")
        flag = "" if ok else "  [unresolved/boundary]"
        plt.title(f"{args.mode}: {sym}={s_c:.3f}+-{dscan:.3f}, |dh|={abs(dh):.1f}, barrier={barrier:.2f} (L only){flag}")
        plt.tight_layout(); plt.savefig(args.plot, dpi=130); print("wrote", args.plot)


if __name__ == "__main__":
    main()
