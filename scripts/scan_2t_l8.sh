#!/bin/bash
# 2T (BT) L=8 (beta,kappa) string-tension grid for the SU(2)->2T phase diagram + triple point.
# Measures, per (beta,kappa) and per hot/cold start: avg_plaquette, acceptance, L_link (Higgs
# order param), fundamental Polyakov, chi(2,2) (the robust small-L confinement proxy = sigma_fund
# upper bound) + its error. Locating the confinement line (chi22 drop) AND the Higgs line (L_link
# rise) over the grid -> the topology + triple point, with the PROPER order parameters, no
# topology assumed. Hot+cold brackets first-order lines.
set -u
cd ~/gh_scan
BIN=build/gh_string
C="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"   # 2T locking quartic couplings (d=7)
MU2=0.113
OUT=scan_2T_L8
mkdir -p "$OUT"
CSV="$OUT/2T_L8.csv"
echo "beta,kappa,start,plaq,acc,Llink,poly,chi22,chi22_err,sigma_line" > "$CSV"
BETAS="0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0"
KAPPAS="0 1 2 3 4 6 8 12"
NTHERM=100; NMEAS=400; NMD=14

run_one() {
  b="$1"; k="$2"; start="$3"
  COLD=""; [ "$start" = cold ] && COLD="GH_COLD=1"
  seed=$(( (RANDOM % 9000) + 100 ))
  out=$(env GH_FROZEN=1 $COLD OMP_NUM_THREADS=3 "$BIN" 6 8 "$b" "$k" "$MU2" "$C" "$NTHERM" "$NMEAS" "$NMD" 1.0 "$seed" 2>/dev/null)
  plaq=$(printf '%s\n' "$out" | grep -m1 'avg_plaquette =' | grep -oE '[-0-9.]+' | head -1)
  acc=$( printf '%s\n' "$out" | grep -m1 'acceptance'      | grep -oE '[0-9.]+'   | head -1)
  Ll=$(  printf '%s\n' "$out" | grep -m1 'L_link'          | grep -oE '[-0-9.]+'  | head -1)
  poly=$(printf '%s\n' "$out" | grep -m1 'polyakov_fund'   | grep -oE '[-0-9.]+'  | head -1)
  c22line=$(printf '%s\n' "$out" | grep -m1 'chi(2,2) =')
  c22v=$(echo "$c22line" | grep -oE '[-0-9.]+ \+/- [-0-9.]+' | head -1 | awk '{print $1}')
  c22e=$(echo "$c22line" | grep -oE '[-0-9.]+ \+/- [-0-9.]+' | head -1 | awk '{print $3}')
  sline=$(printf '%s\n' "$out" | grep -m1 -E 'sigma_fund' | tr ',\n' ';  ')
  echo "$b,$k,$start,${plaq:-NA},${acc:-NA},${Ll:-NA},${poly:-NA},${c22v:-NA},${c22e:-NA},\"${sline:-NA}\"" >> "$CSV"
  echo "done b=$b k=$k $start  chi22=${c22v:-NA} Llink=${Ll:-NA} acc=${acc:-NA}" >> "$OUT/scan.log"
}
export -f run_one
export BIN C MU2 OUT CSV NTHERM NMEAS NMD

echo "START $(date)  $(echo $BETAS | wc -w)x$(echo $KAPPAS | wc -w)x2 = $(( $(echo $BETAS|wc -w)*$(echo $KAPPAS|wc -w)*2 )) runs" > "$OUT/scan.log"
for b in $BETAS; do for k in $KAPPAS; do for s in hot cold; do echo "$b $k $s"; done; done; done \
  | xargs -P 9 -n 3 bash -c 'run_one "$@"' _
echo "ALL DONE $(date)" >> "$OUT/scan.log"
