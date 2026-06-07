#!/bin/bash
# L=12 pure-gauge SU(2) sigma_fund(beta): the CLEAN kappa=0 confinement-line anchor (real Creutz
# plateau over R={2,3,4}, vs the chi22-upper-bound-only at L=8) AND a validation that the L=12
# sigma plateau is clean -- the fundamental Wilson loop is rep-independent, so this de-risks the
# L=12 SU(2)->2T refinement stage. Pure gauge (kappa=0) -> matter decoupled -> cheap; hot start
# (4D SU(2) has no bulk transition). Runs on loranne (4 cores).
set -u
cd ~/gh_scan
BIN=build/gh_string
OUT=scan_puregauge_L12
mkdir -p "$OUT/out"
CSV="$OUT/puregauge_L12.csv"
echo "beta,plaq,acc,chi22,chi33,chi44,sigma_fund_line" > "$CSV"
BETAS="1.0 1.2 1.5 1.8 2.2 2.6"
NTHERM=150; NMEAS=600; NMD=16

run_one() {
  b="$1"
  seed=$(( (RANDOM % 9000) + 100 ))
  o=$(env OMP_NUM_THREADS=2 "$BIN" adj 12 "$b" 0 1 auto "$NTHERM" "$NMEAS" "$NMD" 1.0 "$seed" 2>/dev/null)
  echo "$o" > "$OUT/out/b$b.out"
  plaq=$(echo "$o" | grep -m1 'avg_plaquette =' | grep -oE '[-0-9.]+' | head -1)
  acc=$( echo "$o" | grep -m1 'acceptance'      | grep -oE '[0-9.]+'  | head -1)
  cval() { echo "$o" | grep -m1 "chi($1,$1)" | grep -oE '[-0-9.]+ \+/- [-0-9.]+' | head -1 | awk '{print $1}'; }
  c22=$(cval 2); c33=$(cval 3); c44=$(cval 4)
  sl=$(echo "$o" | grep -m1 -E 'sigma_fund' | tr ',\n' ';  ')
  echo "$b,${plaq:-NA},${acc:-NA},${c22:-NA},${c33:-NA},${c44:-NA},\"${sl:-NA}\"" >> "$CSV"
  echo "done b=$b c22=${c22:-NA} c33=${c33:-NA} c44=${c44:-NA} acc=${acc:-NA}" >> "$OUT/scan.log"
}
export -f run_one; export BIN OUT CSV NTHERM NMEAS NMD
echo "START $(date)  L=12 pure gauge, betas: $BETAS" > "$OUT/scan.log"
for b in $BETAS; do echo "$b"; done | xargs -P 2 -n 1 bash -c 'run_one "$@"' _
echo "ALL DONE $(date)" >> "$OUT/scan.log"
