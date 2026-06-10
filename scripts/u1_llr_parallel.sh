#!/bin/bash
# Level-2 parallel 2D-LLR for the U(1)+charge-q Higgs model: one sequential seed pass
# (drive-in + dump per-cell configs), then NC parallel solve workers over disjoint cell
# ranges, then combine into one reconstruction-ready file. Reproduces `u1_llr` full-mode
# output but parallelizes the (expensive) Robbins-Monro solve OVER CELLS (the constrained
# sweep is inherently serial over sites, so cell-level //ism is the parallel axis).
#
# Build the binary first:   make build/u1_llr      (or `make NDIM=3 build/u1_llr`)
#
# Usage:
#   u1_llr_parallel.sh <outdir> <NC> <L> <q> <lambda> \
#       <Atop> <Abot> <step1> <hw1> <c0> <c1> <c2> <step2> <hw2> <nperp> \
#       <a0> <seed> <K> <NRM>
# Produces <outdir>/combined.out  (feed to scripts/u1_llr_reconstruct.py).
#
# Get the ridge (c0,c1,c2,Atop,Abot) first with:
#   build/u1_llr presample L q lambda npts beta0 beta1 kappa0 kappa1 \
#       | scripts/u1_llr_ridge_presample.py
set -e
OUT=$1; NC=$2; shift 2
L=$1; Q=$2; LAM=$3; shift 3
Atop=$1; Abot=$2; step1=$3; hw1=$4; c0=$5; c1=$6; c2=$7
step2=$8; hw2=$9; nperp=${10}; a0=${11}; seed=${12}; K=${13}; NRM=${14}
LLR=${LLR:-./build/u1_llr}
mkdir -p "$OUT/cfg"

echo "[seed] sequential drive-in -> $OUT/cfg"
"$LLR" seed "$L" "$Q" "$LAM" "$Atop" "$Abot" "$step1" "$hw1" \
    "$c0" "$c1" "$c2" "$step2" "$hw2" "$nperp" "$a0" "$seed" "$OUT/cfg" \
    > "$OUT/seed.log" 2> "$OUT/seed.err"
grep '^CELL:' "$OUT/seed.log" > "$OUT/manifest.txt"
N=$(wc -l < "$OUT/manifest.txt")
echo "[seed] $N cells"

echo "[solve] $NC workers"
per=$(( (N + NC - 1) / NC ))
for w in $(seq 0 $((NC-1))); do
  lo=$((w*per)); hi=$(((w+1)*per))
  [ "$lo" -ge "$N" ] && break
  OMP_NUM_THREADS=1 "$LLR" solve "$L" "$Q" "$LAM" "$OUT/manifest.txt" \
      "$lo" "$hi" "$a0" "$K" "$NRM" "$((seed+1000+w))" \
      > "$OUT/solve_$w.out" 2>> "$OUT/solve.err" &
done
wait

grep '^LLR2D' "$OUT/seed.log" > "$OUT/combined.out"
cat "$OUT"/solve_*.out | grep -E '^(ANE2|RMDIAG):' >> "$OUT/combined.out"
echo "[done] $(grep -c '^ANE2:' "$OUT/combined.out") cells -> $OUT/combined.out"
echo "       reconstruct: scripts/u1_llr_reconstruct.py $OUT/combined.out --q $Q --triple"
