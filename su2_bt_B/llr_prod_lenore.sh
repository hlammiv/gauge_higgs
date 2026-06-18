#!/bin/bash
# SU(2)->2T (BT) LLR production: precise first-order beta_f(kappa) at L=8, kappa=4/5/6/8 (LENORE).
# Two-sided annealing LLR (src/su2_llr.cpp, harmonic restraint, momentum-fixed). Maxwell equal-area on
# a1(e) -> beta_f (scripts/su2_llr_betaf.py). Per-kappa nmd: cold start needs >=20 at k6, more at k8.
# 4 jobs x OMP=7 = 28 threads (fits 32).
cd ~/higgs_gauge
mkdir -p su2_bt_B
FC="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
nmd_for(){ case "$1" in 4) echo 18;; 5) echo 20;; 6) echo 22;; 8) echo 26;; *) echo 20;; esac; }
for k in 4 5 6 8; do
  NMD=$(nmd_for "$k")
  out="su2_bt_B/llr_L8_k${k}.out"
  [ -s "$out" ] && grep -q "^0.4" "$out" && continue
  env GH_FROZEN=1 OMP_NUM_THREADS=7 ./build/su2_llr 6 8 "$k" 0.113 "$FC" \
      0.06 0.48 16 0.7 25 20 0.5 "$NMD" 35 12 $((100+k)) > "$out" 2>&1 &
  sleep 0.5
done
wait
echo "LLR production done: kappa 4,5,6,8 at L=8"
