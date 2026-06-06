#!/bin/bash
# Pure-gauge SU(2) sigma_fund(beta) reference, L=8 (the kappa=0 anchor of the confinement line
# for the SU(2)->2T phase diagram, + a cross-check vs known pure-SU(2) numbers). kappa=0 decouples
# the matter, so this is cheap. Pure SU(2) in 4D is confining at all beta (smooth strong->weak
# crossover, NO bulk transition) -> hot start only, no hysteresis. Runs on loranne (4 cores).
set -u
cd ~/gh_scan
BIN=build/gh_string
OUT=scan_puregauge_L8
mkdir -p "$OUT"
CSV="$OUT/puregauge_L8.csv"
echo "beta,plaq,acc,chi22,chi22_err,sigma_line" > "$CSV"
BETAS="0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0"
NTHERM=120; NMEAS=500; NMD=14

run_one() {
  b="$1"
  seed=$(( (RANDOM % 9000) + 100 ))
  out=$(env OMP_NUM_THREADS=2 "$BIN" adj 8 "$b" 0 1 auto "$NTHERM" "$NMEAS" "$NMD" 1.0 "$seed" 2>/dev/null)
  plaq=$(printf '%s\n' "$out" | grep -m1 'avg_plaquette =' | grep -oE '[-0-9.]+' | head -1)
  acc=$( printf '%s\n' "$out" | grep -m1 'acceptance'      | grep -oE '[0-9.]+'   | head -1)
  c22line=$(printf '%s\n' "$out" | grep -m1 'chi(2,2) =')
  c22v=$(echo "$c22line" | grep -oE '[-0-9.]+ \+/- [-0-9.]+' | head -1 | awk '{print $1}')
  c22e=$(echo "$c22line" | grep -oE '[-0-9.]+ \+/- [-0-9.]+' | head -1 | awk '{print $3}')
  sline=$(printf '%s\n' "$out" | grep -m1 -E 'sigma_fund' | tr ',\n' ';  ')
  echo "$b,${plaq:-NA},${acc:-NA},${c22v:-NA},${c22e:-NA},\"${sline:-NA}\"" >> "$CSV"
  echo "done b=$b chi22=${c22v:-NA} acc=${acc:-NA}" >> "$OUT/scan.log"
}
export -f run_one; export BIN OUT CSV NTHERM NMEAS NMD
echo "START $(date)  $(echo $BETAS | wc -w) betas (pure gauge, hot)" > "$OUT/scan.log"
for b in $BETAS; do echo "$b"; done | xargs -P 2 -n 1 bash -c 'run_one "$@"' _
echo "ALL DONE $(date)" >> "$OUT/scan.log"
