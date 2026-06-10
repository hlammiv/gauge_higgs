#!/usr/bin/env python3
"""1D-A LLR reconstruction + validation at FIXED kappa.

Reads the `ANE1: A0 a1 Bmean sd1 sdB n` lines from `u1_llr slab` (which windows only
the gauge action A and runs the matter physically at fixed kappa, measuring the
physical <B>(A) per cell). Reconstructs ln rho_kappa(A) = integral a1 dA, then
reweights in beta to <A>(beta) and <B>(beta) and compares to the fixed-kappa presample
(direct-HMC ground truth).

Usage:  u1_llr_slab_validate.py <slab.out> <presample.out>
"""
import sys
import numpy as np


def parse_slab(fname):
    A, a1, B = [], [], []
    kappa = None
    for l in open(fname):
        if l.startswith("LLR1DA"):
            for tok in l.replace("|", " ").split():
                if tok.startswith("kappa="):
                    kappa = float(tok.split("=", 1)[1])
        elif l.startswith("ANE1:"):
            t = l.split()
            A.append(float(t[1])); a1.append(float(t[2])); B.append(float(t[3]))
    o = np.argsort(A)
    return np.array(A)[o], np.array(a1)[o], np.array(B)[o], kappa


def softmax(w):
    w = w - np.max(w); e = np.exp(w); return e / e.sum()


def main():
    slab, presamp = sys.argv[1], sys.argv[2]
    A, a1, B, kappa = parse_slab(slab)
    if len(A) < 3:
        sys.exit(f"FAIL: only {len(A)} slab cells")
    # ln rho_kappa(A) = integral a1 dA  (a1 = d ln rho / dA), cumulative trapezoid
    lnrho = np.zeros_like(A)
    lnrho[1:] = np.cumsum(0.5 * (a1[1:] + a1[:-1]) * np.diff(A))
    print(f"# 1D-A slab: {len(A)} cells, kappa={kappa}, A in [{A.min():.1f},{A.max():.1f}], "
          f"a1->beta in [{a1.min():.3f},{a1.max():.3f}]")

    print("\n# beta  (kappa fixed) |   <A>_HMC   <A>_LLR    dA%  |   <B>_HMC   <B>_LLR    dB%")
    dAs, dBs = [], []
    for line in open(presamp):
        if not line.startswith("PRESAMPLE:"):
            continue
        t = line.split(); beta = float(t[1]); kap = float(t[2]); Ah = float(t[3]); Bh = float(t[4])
        if kappa is not None and abs(kap - kappa) > 1e-6:
            continue                      # only the rows at this slab's kappa
        p = softmax(lnrho - beta * A)
        Al = float((A * p).sum()); Bl = float((B * p).sum())
        dA = 100 * (Al - Ah) / Ah if Ah else float("nan")
        dB = 100 * (Bl - Bh) / Bh if Bh else float("nan")
        dAs.append(abs(dA)); dBs.append(abs(dB))
        print(f"{beta:5.2f}              | {Ah:9.1f} {Al:9.1f} {dA:+6.1f}% | {Bh:9.1f} {Bl:9.1f} {dB:+6.1f}%")
    if dAs:
        print(f"\n# GATE: mean|dA| {np.mean(dAs):.1f}%  mean|dB| {np.mean(dBs):.1f}%   "
              f"(both <~5-10% => 1D-A LLR reproduces direct HMC at this kappa)")
    else:
        print("\n# no presample rows matched this kappa")


if __name__ == "__main__":
    main()
