#!/bin/bash
# SU(2)->2T (BT) -- trace the FIRST-ORDER freezing line at high kappa, LENORE.
# The L=8 coarse grid found a chi_plaq spike (47 at k=6,b=2.4; coexistence, tau_int~27 => 1st order). This
# refines beta in {2.2..2.6} at kappa {4,5,6,8} to locate beta_f(kappa) and bracket the critical endpoint
# kappa* (k=4 still crossover, k>=5 1st-order). LONGER runs (meas=600) for the huge autocorr near
# coexistence; per-kappa nmd; center-flip always. NOTE: HMC cannot fully tunnel this line (see memory) ->
# treat as hysteresis-bracketed; precise beta_f is an LLR job. OMP=4, cap 8.
cd ~/higgs_gauge
mkdir -p su2_bt_B
FC="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
CAP=8; i=0
nmd_for(){ case "$1" in 4) echo 12;; 5) echo 14;; 6) echo 16;; 8) echo 20;; *) echo 14;; esac; }
for k in 4 5 6 8; do
  NMD=$(nmd_for "$k")
  for b in 2.2 2.3 2.4 2.5 2.6; do
    out="su2_bt_B/B_L8_b${b}_k${k}.out"      # standard name (plotter reads B_L8_*); 600-meas supersedes coarse
    [ -s "$out" ] && grep -q "600 meas" "$out" && continue
    while [ "$(pgrep -c hmc_higgs)" -ge "$CAP" ]; do sleep 5; done
    seed=$((5000+i)); i=$((i+1))
    env GH_GUPD=2,0,2 GH_FROZEN=1 OMP_NUM_THREADS=4 ./build/hmc_higgs_multi 6 8 "$b" "$k" 0.113 "$FC" 150 600 "$NMD" 1 "$seed" > "$out" 2>&1 &
    sleep 0.3
  done
done
wait
echo "BT first-order-line trace done: $i jobs"
