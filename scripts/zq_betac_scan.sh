#!/bin/bash
# U(1)+charge-q -> Z_q residual: does the large-kappa Higgs route reproduce the KNOWN pure-Z_q
# freezing beta_c at L=2? Abelian analog of the SU(2)->2T bug-hunt. Scan beta through the
# expected Z_q transition(s) at large kappa, hot+cold (bracket first-order), per q.
# For q>=5 expect an intermediate Coulomb phase (TWO transitions): confined / Coulomb / frozen.
# Run from the dir containing build/u1_scan. Set QS (q-list) per machine.
set -u
BIN=build/u1_scan
QS="${QS:?set QS (space-separated q values)}"
KAPPAS="${KAPPAS:-4 8 16}"          # kappa-convergence: does beta_c(kappa) -> known Z_q value?
LAM="${LAM:-1.0}"                    # Higgs quartic (|phi| ~ O(1); large kappa freezes theta -> Z_q)
OUT="${OUT:-zq_betac}"; mkdir -p "$OUT"
BMIN="${BMIN:-0.2}"; BMAX="${BMAX:-2.2}"; NB="${NB:-21}"   # step ~0.1; brackets Z_q beta_c (0.4..1.5) + q>=5 second line
NTH="${NTH:-3000}"; NME="${NME:-2000}"; PAR="${PAR:-8}"

run() {
  q="$1"; k="$2"; st="$3"
  od="$OUT/q${q}_k${k}_${st}"
  E=""; [ "$st" = cold ] && E="U1_COLD=1"
  # args: L bmin bmax nb kmin kmax nk lambda q  ntherm nmeas nmd tau seed measure_every outdir autotune n_scalar
  env $E OMP_NUM_THREADS=2 "$BIN" 2 "$BMIN" "$BMAX" "$NB" "$k" "$k" 1 "$LAM" "$q" \
      "$NTH" "$NME" 20 1.0 12345 1 "$od" 1 1 >/dev/null 2>&1
  echo "done q=$q k=$k $st -> $od/summary.csv ($(tail -n +1 "$od/summary.csv" 2>/dev/null | grep -c '^[0-9.]') rows)"
}
export -f run; export BIN LAM OUT BMIN BMAX NB NTH NME
echo "START $(date)  QS=[$QS] KAPPAS=[$KAPPAS] beta=$BMIN..$BMAX($NB) lambda=$LAM" | tee "$OUT/scan.log"
for q in $QS; do for k in $KAPPAS; do for st in hot cold; do echo "$q $k $st"; done; done; done \
  | xargs -P "$PAR" -n 3 bash -c 'run "$@"' _ | tee -a "$OUT/scan.log"
echo "ALL DONE $(date)" | tee -a "$OUT/scan.log"
