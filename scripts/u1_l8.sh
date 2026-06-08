#!/bin/bash
# U(1)+charge-q L=8 (beta,kappa) cartography (#28 triple point / #26 intermediate-Coulomb).
# Parallelize by beta: each (q,beta) is one u1_scan invocation scanning the kappa-column.
# Run from the dir containing build/u1_scan. Set QS per machine.
set -u
BIN=build/u1_scan
QS="${QS:?set QS (space-separated q values)}"
BETAS="${BETAS:-0.70 0.80 0.90 1.00 1.10 1.20 1.30}"   # brackets U(1) beta_c~1.01
KMIN="${KMIN:-0.0}"; KMAX="${KMAX:-0.6}"; NK="${NK:-7}"   # triple-point kappa window
LAM="${LAM:-1.0}"; L="${L:-8}"
NTH="${NTH:-1000}"; NME="${NME:-3000}"; NMD="${NMD:-10}"; NSC="${NSC:-4}"   # autotune + multi-timescale
OUT="${OUT:-u1_l8}"; mkdir -p "$OUT"; PAR="${PAR:-7}"
run(){ q="$1"; b="$2"; od="$OUT/q${q}_b${b}";
  # args: L bmin bmax nb kmin kmax nk lambda q  ntherm nmeas nmd tau seed measure_every outdir autotune n_scalar
  OMP_NUM_THREADS="${OMPT:-2}" "$BIN" "$L" "$b" "$b" 1 "$KMIN" "$KMAX" "$NK" "$LAM" "$q" \
    "$NTH" "$NME" "$NMD" 1.0 12345 1 "$od" 1 "$NSC" >/dev/null 2>&1
  echo "done q=$q b=$b -> $od/summary.csv"; }
export -f run; export BIN KMIN KMAX NK LAM L NTH NME NMD NSC OUT OMPT
echo "START $(date) QS=[$QS] BETAS=[$BETAS] kappa=$KMIN..$KMAX($NK) L=$L" | tee "$OUT/scan.log"
for q in $QS; do for b in $BETAS; do echo "$q $b"; done; done | xargs -P "$PAR" -n 2 bash -c 'run "$@"' _ | tee -a "$OUT/scan.log"
echo "ALL DONE $(date)" | tee -a "$OUT/scan.log"
