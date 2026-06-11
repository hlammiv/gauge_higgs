#!/usr/bin/env bash
# Emit the frozen-U(1) campaign job queue (one u1_frozen invocation per line), in PRIORITY order:
#   validation -> phase-diagram point grids -> hysteresis (first-order detect) -> muca (tunnel/order).
# Pipe into split + xargs -P. Each line is consumed by scripts/u1f_job.sh.
set -u

# ---- 1. VALIDATION (do first) ----
# pure-gauge beta_c ~ 1.01 (kappa=0): scan beta at several L
for L in 6 8 12; do for b in 0.90 0.95 1.00 1.02 1.05 1.10; do
  echo "point $L $b 0.0 2 5000 2000 7"
done; done
# q=2 deep-kappa -> Z_2 self-dual (beta_eff ~ 0.4407): scan beta at large kappa
for k in 4.0 8.0; do for b in 0.20 0.30 0.40 0.44 0.45 0.50 0.60 0.70; do
  echo "point 8 $b $k 2 5000 2000 7"
done; done

# ---- 2. PHASE-DIAGRAM POINT GRID (observables A,B,plaq,<cos> over the plane) ----
BETAS="0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20"
KAPS="0.10 0.20 0.30 0.50 0.70 1.00 1.50 2.00 3.00"
for q in 2 3 4 5 6 8; do for b in $BETAS; do for k in $KAPS; do
  echo "point 8 $b $k $q 4000 1500 7"
done; done; done
# FSS hint: q=2,4 at L=12 on a coarser grid
for q in 2 4; do for b in 0.50 0.70 0.90 1.10; do for k in 0.20 0.50 1.00 2.00; do
  echo "point 12 $b $k $q 4000 1500 7"
done; done; done

# ---- 3. HYSTERESIS (Higgs line; gap => first order) ----
for q in 2 3 4 5 6 8; do for b in 0.70 0.90; do
  echo "hyst 8 $b $q 0.10 2.00 20 600 1500 7 3"
done; done

# ---- 4. MULTICANONICAL (tunnel/order) -- slower, run last ----
# muca-in-B (deep-Higgs matter jump): B in [0, 2*vol*D] = [0,32768] at L=8
for q in 2 4 6 8; do
  echo "muca 8 0.90 0.50 $q 0 32768 64 150000 2000 7 3 u1f_campaign/mucaB_q${q}"
done
# muca-in-A (large-beta Z_q ordering barrier): A in [0, 2*n_plaq] = [0,49152]; deep Higgs kappa=2
for q in 2 4 6 8; do
  echo "mucaA 8 0.90 2.0 $q 0 24576 48 150000 2000 7 3 u1f_campaign/mucaA_q${q}"
done
