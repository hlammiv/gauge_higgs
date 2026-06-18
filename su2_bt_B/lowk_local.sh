#!/bin/bash
# SU(2)->2T (BT) Stage B EXTENSION -- LOW-kappa SCALING REGIME, L=8, LOCAL.
# Fills kappa in {0,0.1,0.2,0.3,0.4} (between 0 and the existing 0.5) with FINE beta {1.6..2.6 step 0.1}
# and BETTER stats (therm 100, meas 400) to sharpen the freezing line beta_f(kappa) where it emerges from
# the pure-SU(2) bulk crossover (the "scaling regime at low kappa"). Closed-form D active; OMP=1 x many.
cd ~/Desktop/QC/higgs+gauge
mkdir -p su2_bt_B
FC="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
CAP=16; i=0
for k in 0 0.1 0.2 0.3 0.4; do
  for b in 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6; do
    out="su2_bt_B/B_L8_b${b}_k${k}.out"
    # overwrite the coarse/short k=0 points; keep finished fine ones
    [ -s "$out" ] && grep -q acceptance "$out" && grep -q "n_meas=400\|400 meas" "$out" && continue
    while [ "$(pgrep -c hmc_higgs)" -ge "$CAP" ]; do sleep 5; done
    seed=$((3000+i)); i=$((i+1))
    env GH_FROZEN=1 OMP_NUM_THREADS=1 ./build/hmc_higgs_multi 6 8 "$b" "$k" 0.113 "$FC" 100 400 8 1 "$seed" > "$out" 2>&1 &
    sleep 0.2
  done
done
wait
echo "BT low-kappa scaling-regime grid done: $i jobs"
