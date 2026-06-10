#!/bin/bash
# Parallel RECTANGULAR 2D-LLR for the first-order-gap region. One sequential seed pass
# (two-pass meet-in-the-gap raster -> drive each cell in -> dump config + manifest), then
# NC parallel solve workers do the (expensive) K*NRM Robbins-Monro over disjoint cell
# ranges. Makes K=100/NRM=200 tractable (the rect raster is sequential; the RM is the cost).
#
# Set BOTH ridges first (a1~beta(A), a2~kappa(A)) from the grid presample so a2 starts near
# kappa (else a2 stays pinned at its init -> dead kappa-axis):
#   export U1_LLR_BRIDGE="d0 d1 d2"   U1_LLR_KRIDGE="e0 e1 e2"
#
# Usage:
#   u1_llr_rect_parallel.sh <outdir> <NC> <L> <q> <lambda> \
#       <Atop> <Abot> <step1> <hw1> <Bmin> <Bmax> <step2> <hw2> <a0> <seed> <K> <NRM>
# Produces <outdir>/combined.out  (feed to scripts/u1_llr_rect_reconstruct.py).
set -e
OUT=$1; NC=$2; shift 2
L=$1; Q=$2; LAM=$3; shift 3
Atop=$1; Abot=$2; step1=$3; hw1=$4; Bmin=$5; Bmax=$6; step2=$7; hw2=$8
a0=$9; seed=${10}; K=${11}; NRM=${12}
LLR=${LLR:-./build/u1_llr}
mkdir -p "$OUT/cfg"

echo "[rectseed] two-pass meet-in-the-gap drive-in -> $OUT/cfg"
"$LLR" rectseed "$L" "$Q" "$LAM" "$Atop" "$Abot" "$step1" "$hw1" "$Bmin" "$Bmax" "$step2" "$hw2" \
    "$a0" "$seed" "$OUT/cfg" > "$OUT/seed.log" 2> "$OUT/seed.err"
grep '^CELL:' "$OUT/seed.log" > "$OUT/manifest.txt"
N=$(wc -l < "$OUT/manifest.txt")
echo "[rectseed] $N cells dumped"

echo "[solve] $NC workers (K=$K NRM=$NRM), env ridges: BRIDGE=[${U1_LLR_BRIDGE:-unset}] KRIDGE=[${U1_LLR_KRIDGE:-unset}]"
per=$(( (N + NC - 1) / NC ))
for w in $(seq 0 $((NC-1))); do
  lo=$((w*per)); hi=$(((w+1)*per))
  [ "$lo" -ge "$N" ] && break
  OMP_NUM_THREADS=1 "$LLR" solve "$L" "$Q" "$LAM" "$OUT/manifest.txt" \
      "$lo" "$hi" "$a0" "$K" "$NRM" "$((seed+1000+w))" \
      > "$OUT/solve_$w.out" 2>> "$OUT/solve.err" &
done
wait

grep '^LLR2DRECT' "$OUT/seed.log" > "$OUT/combined.out"
grep '^NPLAQ:\|^NLINKS:' "$OUT/seed.log" >> "$OUT/combined.out"
cat "$OUT"/solve_*.out | grep -E '^(ANE2|RMDIAG):' >> "$OUT/combined.out"
echo "[done] $(grep -c '^ANE2:' "$OUT/combined.out") cells -> $OUT/combined.out"
echo "       reconstruct: scripts/u1_llr_rect_reconstruct.py $OUT/combined.out --q $Q"
