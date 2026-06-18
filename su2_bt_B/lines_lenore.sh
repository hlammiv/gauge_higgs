#!/bin/bash
# SU(2)->2T (BT) Stage B, L=12 FREEZING-LINE -- runs LENORE (heavier; OMP=4, closed-form D).
# chi_plaq FSS near the freezing transition (beta~1.8-2.6) at kappa=0.5/1.0/2.0; locates beta_f(kappa)
# and FSS-confirms it vs the L=8 grid. nmd=16 (cal: L=12 k=2 nmd=16 -> acc 0.76). Center-flip at k>=1.
# Cap 8 concurrent x OMP=4 = 32 threads (full box, no oversubscribe).
cd ~/higgs_gauge
mkdir -p su2_bt_B
FC="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
CAP=8; i=0
for k in 0.5 1.0 2.0; do
  G=""; awk "BEGIN{exit !($k>=1)}" && G='GH_GUPD=2,0,2'
  for b in 1.8 2.0 2.1 2.2 2.3 2.4 2.6; do
    out="su2_bt_B/B_L12_b${b}_k${k}.out"
    [ -s "$out" ] && grep -q acceptance "$out" && continue
    while [ "$(pgrep -c hmc_higgs)" -ge "$CAP" ]; do sleep 5; done
    seed=$((2000+i)); i=$((i+1))
    env $G GH_FROZEN=1 OMP_NUM_THREADS=4 ./build/hmc_higgs_multi 6 12 "$b" "$k" 0.113 "$FC" 75 225 16 1 "$seed" > "$out" 2>&1 &
    sleep 0.3
  done
done
wait
echo "BT Stage-B L=12 freezing-line done: $i jobs"
