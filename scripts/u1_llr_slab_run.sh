#!/bin/bash
# 1D-A LLR validation at FIXED kappa: window only A, matter physical (heavy phi sampling).
# fixed-kappa presample (ground truth, vary beta) -> beta(A) ridge -> slab tile A -> reweight gate.
# Usage: u1_llr_slab_run.sh <L> <kappa> [q NA K NRM]
set -e
L="${1:?L}"; KAPPA="${2:?kappa}"; Q="${3:-2}"; NA="${4:-15}"; K="${5:-60}"; NRM="${6:-120}"; LAM=0.5
BIN=./build/u1_llr; OUT="llr_slab_L${L}_q${Q}_k${KAPPA}"; mkdir -p "$OUT"
# heavy matter sampling (validated: n_site=8 tau_s=1.0 n_md_s=20 -> <B> within ~3%)
DELTA0=0.6; NHITG=5; NOVER=2; TAUS=1.0; NMDS=20; NSITE=8; SITEW=0.3

echo "[$(date +%T)] L=$L q=$Q kappa=$KAPPA  fixed-kappa presample (vary beta)..."
OMP_NUM_THREADS=8 "$BIN" presample "$L" "$Q" "$LAM" 8 0.80 1.50 "$KAPPA" "$KAPPA" 150 400 15 1.0 7 \
    > "$OUT/presamp.out" 2>/dev/null
python3 scripts/u1_llr_ridge_presample.py "$OUT/presamp.out" > "$OUT/ridge.out" 2>/dev/null
eval "$(grep '^export' "$OUT/ridge.out")"; export U1_LLR_BRIDGE         # a1~beta(A) seed
Atop=$(grep -m1 -oP 'Atop=\K[-0-9.eE+]+' "$OUT/ridge.out"); Abot=$(grep -m1 -oP 'Abot=\K[-0-9.eE+]+' "$OUT/ridge.out")
step1=$(awk -v at="$Atop" -v ab="$Abot" -v n="$NA" 'BEGIN{printf "%.4f", (at-ab)/(n-1)}')
hw1="$step1"
echo "[$(date +%T)] A=[$Abot,$Atop] step1=$step1 hw1=$hw1 BRIDGE=$U1_LLR_BRIDGE"
echo "[$(date +%T)] slab tile (NA=$NA K=$K NRM=$NRM heavy-matter n_site=$NSITE tau_s=$TAUS n_md_s=$NMDS)..."
OMP_NUM_THREADS=8 "$BIN" slab "$L" "$Q" "$LAM" "$KAPPA" "$Atop" "$Abot" "$step1" "$hw1" 1.0 7 "$K" "$NRM" \
    "$DELTA0" "$NHITG" "$NOVER" "$TAUS" "$NMDS" "$NSITE" "$SITEW" > "$OUT/slab.out" 2>"$OUT/slab.err"
echo "[$(date +%T)] validating vs HMC ground truth..."
python3 scripts/u1_llr_slab_validate.py "$OUT/slab.out" "$OUT/presamp.out" | tee "$OUT/validate.out"
echo "[$(date +%T)] DONE L=$L kappa=$KAPPA -> $OUT/validate.out"
