#!/usr/bin/env python3
"""
Ridge-fit helper for the 2D-LLR U(1)+Higgs run. Reads the `PRESAMPLE:` lines emitted by
`u1_llr presample` (plain U1HMC at several (beta,kappa) along the expected ribbon) and
least-squares-fits the curved ridge

    E2(A) = -B(A) = c0 + c1*A + c2*A^2

then prints the (c0,c1,c2) and the [Atop, Abot] band endpoints (with a 10% margin) ready to
paste into the `u1_llr` full / seed CLI. The band runs from the disordered/confined corner
(LARGE A, small B => E2 near 0) to the ordered Higgs/Coulomb corner (SMALL A, large B =>
E2 large-negative); c2 captures the curvature.

Usage:
  build/u1_llr presample L q lambda npts beta0 beta1 kappa0 kappa1 ... | \
      scripts/u1_llr_ridge_presample.py [--margin 0.10] [--plot ridge.png]
or
  scripts/u1_llr_ridge_presample.py presample.out [--margin 0.10]
"""
import sys, argparse
import numpy as np


def parse(stream):
    beta, kappa, A, B, E2 = [], [], [], [], []
    for line in stream:
        if not line.startswith("PRESAMPLE:"):
            continue
        t = line.split()
        # PRESAMPLE: beta kappa A_mean B_mean E2_mean
        beta.append(float(t[1])); kappa.append(float(t[2]))
        A.append(float(t[3])); B.append(float(t[4])); E2.append(float(t[5]))
    return (np.array(beta), np.array(kappa), np.array(A), np.array(B), np.array(E2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="?", default="-")
    ap.add_argument("--margin", type=float, default=0.10,
                    help="fractional A-margin added beyond the sampled [minA,maxA] (default 0.10)")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    stream = sys.stdin if args.file == "-" else open(args.file)
    beta, kappa, A, B, E2 = parse(stream)
    if len(A) < 3:
        sys.exit("need >=3 PRESAMPLE points to fit a quadratic ridge")

    # least-squares quadratic E2 = c0 + c1*A + c2*A^2  (energy ridge)
    c2, c1, c0 = np.polyfit(A, E2, 2)
    resid = E2 - (c0 + c1 * A + c2 * A * A)
    rms = float(np.sqrt(np.mean(resid ** 2)))

    # RM-init ridges so a1~beta(A) and a2~kappa(A) start near the answer on BOTH axes:
    #   beta(A)  = d0 + d1*A + d2*A^2     (-> a1 seed,  env U1_LLR_BRIDGE)
    #   kappa(A) = e0 + e1*A + e2*A^2     (-> a2 seed,  env U1_LLR_KRIDGE)
    # (use a line if only ~3 points; quadratic once there are >=4)
    deg = 2 if len(A) >= 4 else 1
    bcoef = np.polyfit(A, beta, deg)          # numpy returns coeffs HIGH power -> low
    kcoef = np.polyfit(A, kappa, deg)
    # rebuild as low->high [const, lin, quad] padded to 3 for the driver env
    bc = np.zeros(3); bc[:deg + 1] = bcoef[::-1]; d0b, d1b, d2b = bc
    kc = np.zeros(3); kc[:deg + 1] = kcoef[::-1]; e0k, e1k, e2k = kc

    Amin, Amax = float(A.min()), float(A.max())
    span = Amax - Amin
    Atop = Amax + args.margin * span        # disordered / confined end (large A)
    Abot = max(0.0, Amin - args.margin * span)  # ordered end (small A); A>=0

    print(f"# ridge fit over {len(A)} presample points;  E2 = c0 + c1*A + c2*A^2   (RMS resid {rms:.3g})")
    print(f"# A sampled [{Amin:.2f}, {Amax:.2f}]   E2 sampled [{E2.min():.2f}, {E2.max():.2f}]")
    print(f"# beta sampled [{beta.min():.3f}, {beta.max():.3f}]   kappa sampled [{kappa.min():.3f}, {kappa.max():.3f}]  (fit deg {deg})")
    print(f"c0={c0:.6g} c1={c1:.6g} c2={c2:.6g}")
    print(f"Atop={Atop:.4f} Abot={Abot:.4f}")
    print()
    print("# RM-init ridges -- export these so the driver seeds a1~beta(A), a2~kappa(A) per cell:")
    print(f'export U1_LLR_BRIDGE="{d0b:.8g} {d1b:.8g} {d2b:.8g}"')
    print(f'export U1_LLR_KRIDGE="{e0k:.8g} {e1k:.8g} {e2k:.8g}"')
    print()
    print("# paste-ready geometry tail (pick step1/hw1/step2/hw2/nperp per the contract):")
    print(f"#   ... Atop={Atop:.2f} Abot={Abot:.2f} step1=<(Atop-Abot)/(N1-1)> hw1=<step1> "
          f"c0={c0:.5g} c1={c1:.5g} c2={c2:.5g} step2=<hw2> hw2=<max(2*step2,3*sigma_dB)> nperp=3 ...")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = np.linspace(Abot, Atop, 200)
        plt.figure(figsize=(6.5, 5))
        plt.scatter(A, E2, c="k", s=40, label="presample (A, -B)")
        plt.plot(xs, c0 + c1 * xs + c2 * xs * xs, "r-", lw=1.8, label="quadratic ridge")
        plt.axvline(Atop, color="gray", ls=":"); plt.axvline(Abot, color="gray", ls=":")
        plt.xlabel("A = sum_plaq(1-cos)"); plt.ylabel("E2 = -B")
        plt.title("2D-LLR ridge fit (U(1)+Higgs)")
        plt.legend(); plt.tight_layout(); plt.savefig(args.plot, dpi=120)
        print("plot ->", args.plot)


if __name__ == "__main__":
    main()
