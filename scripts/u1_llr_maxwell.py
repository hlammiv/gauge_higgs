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
    """Return (i_lo, i_hi, i_valley) for the two highest separated local maxima of P, or None."""
    pk = [i for i in range(len(P)) if (i == 0 or P[i] >= P[i-1]) and (i == len(P)-1 or P[i] >= P[i+1])]
    pk = [i for i in pk if np.isfinite(P[i])]
    if len(pk) < 2:
        return None
    # take the two highest peaks that are separated by a valley
    pk.sort(key=lambda i: -P[i])
    i1, i2 = sorted(pk[:2])
    iv = i1 + int(np.argmin(P[i1:i2+1]))
    if iv == i1 or iv == i2:
        return None
    return i1, i2, iv


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
    best = None
    rows = []
    for s in scans:
        P = logsumexp_bins(expo(s), Ok, uO)
        tp = two_peaks(P)
        if tp is None:
            rows.append((s, None)); continue
        i1, i2, iv = tp
        dh = P[i1] - P[i2]                     # peak-height difference (0 => equal weight)
        barrier = min(P[i1], P[i2]) - P[iv]    # Delta F (in -ln P units)
        latent = abs(uO[i2] - uO[i1])
        rows.append((s, (dh, barrier, latent, uO[i1], uO[i2])))
        if best is None or abs(dh) < abs(best[1][0]):
            best = (s, (dh, barrier, latent, uO[i1], uO[i2]), P)

    print(f"# {args.mode} transition, fixed {fixed}, scan {args.scan}, {len(cells)} cells, V={V}")
    ndp = sum(1 for s, r in rows if r is not None)
    print(f"# double-peak resolved in {ndp}/{len(scans)} scan points "
          + ("(NONE -> not resolved as first-order at this L/resolution => crossover-like)" if ndp == 0 else ""))
    if best is None:
        print("# NO equal-weight double peak found -> transition NOT numerically resolved here.")
        return
    s_c, (dh, barrier, latent, olo, ohi), P = best
    sym = "kappa_c" if args.mode == "higgs" else "beta_c"
    print(f"# EQUAL-WEIGHT {sym} = {s_c:.4f}  (peak-height diff {dh:+.3f})")
    print(f"# latent heat  Delta{'B' if args.mode=='higgs' else 'A'} = {ohi-olo:.1f}  (peaks at {olo:.0f}, {ohi:.0f})")
    print(f"# free-energy BARRIER Delta F = {barrier:.3f}  (per-volume {barrier/V:.5f}); grows with V => FIRST ORDER")

    if args.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        plt.figure(figsize=(6.5, 4.5))
        Pn = P - P.max()
        plt.plot(uO, Pn, "o-", ms=3)
        plt.axvline(olo, color="g", ls=":"); plt.axvline(ohi, color="r", ls=":")
        plt.xlabel("B (matter)" if args.mode == "higgs" else "A (gauge)")
        plt.ylabel("ln P (reweighted)")
        plt.title(f"{args.mode}: equal-weight {sym}={s_c:.3f}, latent={ohi-olo:.0f}, barrier={barrier:.2f}")
        plt.tight_layout(); plt.savefig(args.plot, dpi=130); print("wrote", args.plot)


if __name__ == "__main__":
    main()
