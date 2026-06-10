#!/usr/bin/env python3
"""Validation gate for a 2D-LLR run: reweight the reconstructed rho(A,B) to the
ribbon (beta,kappa) points and compare <A>,<B> against the direct-HMC presample
(the ground truth), plus a per-cell RM-convergence summary.

Usage:  u1_llr_validate.py <combined.out> <presample.out> [amax]
"""
import sys, os, importlib.util
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("rec", os.path.join(_here, "u1_llr_reconstruct.py"))
rec = importlib.util.module_from_spec(spec); spec.loader.exec_module(rec)


def softmax(w):
    w = w - np.max(w); e = np.exp(w); return e / e.sum()


def main():
    combined, presamp = sys.argv[1], sys.argv[2]
    amax = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0
    cells, geo, sd, rmdiag = rec.parse(combined, amax)
    if len(cells) < 4:
        sys.exit(f"FAIL: only {len(cells)} converged cells (|a|<{amax}) -- tile did not converge")
    idx = rec.grid_index(cells, geo)
    lnrho = rec.integrate(cells, idx, sd)
    A = cells[:, 0]; E2 = cells[:, 1]; B = -E2

    # --- convergence summary from RMDIAG ---
    hw1 = geo.get("hw1", float("nan")); hw2 = geo.get("hw2", float("nan"))
    if rmdiag:
        rA = np.array([d["res_A"] for d in rmdiag.values()])
        rE = np.array([d["res_E2"] for d in rmdiag.values()])
        ah = np.array([d["acc_hmc"] for d in rmdiag.values()])
        nw = np.array([d["nwiden"] for d in rmdiag.values()])
        print(f"# {len(cells)} ANE2 cells kept; {len(rmdiag)} RMDIAG")
        print(f"# CONVERGENCE: res_A med {np.median(rA):.2f} (hw1={hw1:.1f}) frac<0.3hw1 {np.mean(rA<0.3*hw1):.2f} | "
              f"res_E2 med {np.median(rE):.2f} (hw2={hw2:.1f}) frac<0.3hw2 {np.mean(rE<0.3*hw2):.2f}")
        print(f"# acc_hmc med {np.median(ah):.2f} (want ~0.5-0.9) | nwiden==0 frac {np.mean(nw==0):.2f}")
    print(f"# a1->beta range [{cells[:,2].min():.3f},{cells[:,2].max():.3f}] | "
          f"a2->kappa range [{cells[:,3].min():.3f},{cells[:,3].max():.3f}]")

    # --- reweight gate vs presample ground truth ---
    # Use the CONTINUOUS-tile moments (cell_moments): each cell is a finite box [E1±hw1]x[E2±hw2]
    # with the (a1,a2) tilt, so the within-cell mean offset (ub) is included. The crude cell-CENTRE
    # delta drops ub, which for wide E2 cells (hw2~tens) badly biases <B>. <A>=S1, <B>=-S2.
    hw1 = float(geo["hw1"]); hw2 = float(geo["hw2"])
    print("\n# beta  kappa |   <A>_HMC   <A>_LLR    dA%  |   <B>_HMC   <B>_LLR    dB%   (continuous-tile moments)")
    dAs, dBs = [], []
    for line in open(presamp):
        if not line.startswith("PRESAMPLE:"):
            continue
        t = line.split(); beta = float(t[1]); kappa = float(t[2]); Ah = float(t[3]); Bh = float(t[4])
        S1, _, S2, _ = rec.cell_moments(cells, lnrho, hw1, hw2, beta, kappa)
        Al = float(S1); Bl = float(-S2)
        dA = 100 * (Al - Ah) / Ah if Ah else float("nan")
        dB = 100 * (Bl - Bh) / Bh if Bh else float("nan")
        dAs.append(abs(dA)); dBs.append(abs(dB))
        print(f"{beta:5.2f} {kappa:5.3f} | {Ah:9.1f} {Al:9.1f} {dA:+6.1f}% | {Bh:9.1f} {Bl:9.1f} {dB:+6.1f}%")
    print(f"\n# GATE: mean|dA| {np.mean(dAs):.1f}%  mean|dB| {np.mean(dBs):.1f}%   "
          f"(<~5% A and <~10% B => LLR reproduces direct HMC)")


if __name__ == "__main__":
    main()
