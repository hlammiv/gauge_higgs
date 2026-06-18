#!/bin/bash
# SU(2)->2T (BT) -- COLD-start companion for HYSTERESIS bracketing of the 1st-order freezing line, LENORE.
# The hot-start foline runs scatter (HMC can't tunnel the strong 1st-order line -> configs randomly freeze).
# Driver: GH_COLD = ordered start (identity links + aligned phi); default = hot. Where plaquette(cold) and
# plaquette(hot) DIFFER at the same (beta,kappa) = the hysteresis loop -> beta_f is bracketed inside it
# (Bowler PLB104B'81 / Damgaard-Heller NPB324'89 method). Pair with the existing hot B_L8_b*_k{4,5,6,8}.
# meta: writes B_L8cold_b*_k*.out (separate from the hot set). per-kappa nmd; OMP=4, cap 8.
cd ~/higgs_gauge
mkdir -p su2_bt_B
FC="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
CAP=8; i=0
nmd_for(){ case "$1" in 4) echo 12;; 5) echo 14;; 6) echo 16;; 8) echo 20;; *) echo 14;; esac; }
for k in 4 5 6 8; do
  NMD=$(nmd_for "$k")
  for b in 2.1 2.2 2.3 2.4 2.5 2.6; do
    out="su2_bt_B/B_L8cold_b${b}_k${k}.out"
    [ -s "$out" ] && grep -q acceptance "$out" && continue
    while [ "$(pgrep -c hmc_higgs)" -ge "$CAP" ]; do sleep 5; done
    seed=$((6000+i)); i=$((i+1))
    env GH_COLD=1 GH_GUPD=2,0,2 GH_FROZEN=1 OMP_NUM_THREADS=4 ./build/hmc_higgs_multi 6 8 "$b" "$k" 0.113 "$FC" 150 600 "$NMD" 1 "$seed" > "$out" 2>&1 &
    sleep 0.3
  done
done
wait
echo "BT cold-start hysteresis scan done: $i jobs"
