#!/bin/bash
# End-to-end 2D-LLR validation gate for U(1)+charge-q on one lattice size.
# presample (= ground truth) -> ridge fit -> converged parallel tile -> reweight-vs-HMC.
# Usage: u1_llr_validate_run.sh <L> [q] [NC] [K] [NRM]
set -e
L="${1:?L}"; Q="${2:-2}"; NC="${3:-28}"; K="${4:-100}"; NRM="${5:-200}"; LAM=0.5
BIN=./build/u1_llr; OUT="llr_val_L${L}_q${Q}"; mkdir -p "$OUT"
if [ -s "$OUT/presamp.out" ] && [ "$(grep -c '^PRESAMPLE:' "$OUT/presamp.out")" -ge 4 ]; then
  echo "[$(date +%T)] L=$L q=$Q  reusing existing presample ($(grep -c '^PRESAMPLE:' "$OUT/presamp.out") pts)"
else
  echo "[$(date +%T)] L=$L q=$Q  presample (ground truth)..."
  OMP_NUM_THREADS=8 "$BIN" presample "$L" "$Q" "$LAM" 8 0.85 1.45 0.07 0.18 150 400 15 1.0 7 \
      > "$OUT/presamp.out" 2>/dev/null
fi
python3 scripts/u1_llr_ridge_presample.py "$OUT/presamp.out" > "$OUT/ridge.out" 2>/dev/null
eval "$(grep '^export' "$OUT/ridge.out")"            # -> U1_LLR_BRIDGE, U1_LLR_KRIDGE
export U1_LLR_BRIDGE U1_LLR_KRIDGE
# -m1: take ONLY the real fit line, not the paste-ready comment that repeats c0=/Atop=
c0=$(grep -m1 -oP 'c0=\K[-0-9.eE+]+' "$OUT/ridge.out"); c1=$(grep -m1 -oP ' c1=\K[-0-9.eE+]+' "$OUT/ridge.out")
c2=$(grep -m1 -oP ' c2=\K[-0-9.eE+]+' "$OUT/ridge.out")
Atop=$(grep -m1 -oP 'Atop=\K[-0-9.eE+]+' "$OUT/ridge.out"); Abot=$(grep -m1 -oP 'Abot=\K[-0-9.eE+]+' "$OUT/ridge.out")
# geometry from the A and E2 spans: N1=20 ridge cells, band hw2=0.30*E2span, nperp=2
read step1 hw1 step2 hw2 < <(awk -v at="$Atop" -v ab="$Abot" '
  /PRESAMPLE:/{e=$6; if(n++==0){mn=e;mx=e} if(e<mn)mn=e; if(e>mx)mx=e}
  END{ s1=(at-ab)/19.0; span=mx-mn; if(span<1)span=1; h2=0.30*span; if(h2<40)h2=40;
       printf "%.4f %.4f %.4f %.4f\n", s1, s1, h2/2.0, h2 }' "$OUT/presamp.out")
echo "[$(date +%T)] ridge c0=$c0 c1=$c1 c2=$c2  A=[$Abot,$Atop] step1=$step1 hw1=$hw1 step2=$step2 hw2=$hw2"
echo "             BRIDGE=$U1_LLR_BRIDGE  KRIDGE=$U1_LLR_KRIDGE"
echo "[$(date +%T)] tiling (NC=$NC K=$K NRM=$NRM nperp=2)..."
scripts/u1_llr_parallel.sh "$OUT/tile" "$NC" "$L" "$Q" "$LAM" \
    "$Atop" "$Abot" "$step1" "$hw1" "$c0" "$c1" "$c2" "$step2" "$hw2" 2 0.5 7 "$K" "$NRM"
echo "[$(date +%T)] validating vs HMC ground truth..."
python3 scripts/u1_llr_validate.py "$OUT/tile/combined.out" "$OUT/presamp.out" | tee "$OUT/validate.out"
echo "[$(date +%T)] DONE L=$L -> $OUT/validate.out"
