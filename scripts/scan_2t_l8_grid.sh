#!/bin/bash
# SU(2)->2T L=8 (beta,kappa) PHASE-DIAGRAM campaign (publication-grade #23 first pass).
# Per node: sigma_fund proxy chi(2,2)/chi(3,3) (Wilson-loop string tension, the rigorous
# confinement order param at L=8, Rmax=4) + L_link (Higgs) + plaq + fundamental Polyakov.
# hot+cold per node -> brackets the metastable / first-order-like deep-Higgs band.
# Parameterized: set BETAS (space list) per machine; KAPPAS/stats overridable. Run from the
# dir containing build/gh_string (lenore/loranne: ~/gh_scan ; lucia: the repo).
set -u
BIN=build/gh_string
C="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"   # 2T-locking couplings (d=7)
MU2=0.113; L=8
OUT="${OUT:-scan_2t_l8grid}"; mkdir -p "$OUT"
CSV="$OUT/2t_l8_${TAG:-part}.csv"
echo "beta,kappa,start,plaq,acc,Llink,chi_link,poly,chi22,chi33,sigma_line" > "$CSV"
BETAS="${BETAS:?set BETAS (space-separated)}"
KAPPAS="${KAPPAS:-0 1 2 3 4 6 8}"
NTHERM="${NTHERM:-120}"; NMEAS="${NMEAS:-500}"; NMD="${NMD:-16}"; PAR="${PAR:-10}"; OMPT="${OMPT:-3}"

run_one() {
  b="$1"; k="$2"; st="$3"
  E="GH_FROZEN=1"; [ "$st" = cold ] && E="$E GH_COLD=1"
  o=$(env $E OMP_NUM_THREADS=$OMPT "$BIN" 6 "$L" "$b" "$k" "$MU2" "$C" "$NTHERM" "$NMEAS" "$NMD" 1.0 12345 2>/dev/null)
  g(){ echo "$o"|grep -m1 "$1"|grep -oE '[-0-9.]+'|head -1; }
  cv(){ echo "$o"|grep -m1 "^chi($1,$1) ="|grep -oE '[-0-9.]+ \+/- [-0-9.]+'|head -1|awk '{print $1}'; }
  plaq=$(g 'avg_plaquette ='); acc=$(echo "$o"|grep -m1 acceptance|grep -oE '[0-9.]+'|head -1)
  ll=$(g '^L_link'); poly=$(g 'polyakov_fund'); c22=$(cv 2); c33=$(cv 3); cl=$(g 'chi_link =')
  sl=$(echo "$o"|grep -m1 -E 'sigma_fund'|tr ',\n' ';  ')
  echo "$b,$k,$st,${plaq:-NA},${acc:-NA},${ll:-NA},${cl:-NA},${poly:-NA},${c22:-NA},${c33:-NA},\"${sl:-NA}\"" >> "$CSV"
  echo "done b=$b k=$k $st chi22=${c22:-NA} Ll=${ll:-NA} acc=${acc:-NA}" >> "$OUT/scan_${TAG:-part}.log"
}
export -f run_one; export BIN C MU2 L OUT CSV NTHERM NMEAS NMD OMPT KAPPAS
echo "START $(date)  BETAS=[$BETAS] KAPPAS=[$KAPPAS] NMEAS=$NMEAS NMD=$NMD PAR=$PAR" > "$OUT/scan_${TAG:-part}.log"
for b in $BETAS; do for k in $KAPPAS; do for st in hot cold; do echo "$b $k $st"; done; done; done \
  | xargs -P "$PAR" -n 3 bash -c 'run_one "$@"' _
echo "ALL DONE $(date)" >> "$OUT/scan_${TAG:-part}.log"
