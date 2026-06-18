#!/bin/bash
# SU(2)->2T (BT) Stage B, L=8 COARSE (beta,kappa) MAP -- runs LOCAL (light; closed-form D + nmd=8).
# Gauge freezing line via chi_plaq + matter line via chi_link. kappa<=2.5 only (deep-kappa needs high nmd,
# deferred). Center-flip GH_GUPD at kappa>=1. nmd=8 (pilot: nmd=12 -> accept 0.9 at L=8; 8 -> ~0.7-0.8).
# OMP=1 x many (throughput on the 20-core box). Closed-form D auto-active for spin-3.
cd ~/Desktop/QC/higgs+gauge
mkdir -p su2_bt_B
FC="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
CAP=18; i=0
for k in 0 0.5 1.0 1.5 2.0 2.5; do
  GUPD=""; awk "BEGIN{exit !($k>=1)}" && GUPD='GH_GUPD=2,0,2'
  for b in 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8; do
    out="su2_bt_B/B_L8_b${b}_k${k}.out"
    [ -s "$out" ] && grep -q acceptance "$out" && continue
    while [ "$(pgrep -c hmc_higgs)" -ge "$CAP" ]; do sleep 5; done
    seed=$((1000+i)); i=$((i+1))
    env $GUPD GH_FROZEN=1 OMP_NUM_THREADS=1 ./build/hmc_higgs_multi 6 8 "$b" "$k" 0.113 "$FC" 50 150 8 1 "$seed" > "$out" 2>&1 &
    sleep 0.2
  done
done
wait
echo "BT Stage-B L=8 grid done: $i jobs"
