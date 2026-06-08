#!/bin/bash
# SU(2)->2T (beta,kappa) PHASE-DIAGRAM grid, L=2, metastability-aware: hot AND cold start per node.
# Observables: <plaq> (gauge axis: confined low / weak-coupling high) + L_link (matter axis: Higgs
# when high). hot-vs-cold gap flags the metastable / first-order-like region (deep Higgs). frozen 2T.
set -u
cd ~/gh_scan
BIN=build/gh_string
C="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"   # 2T-locking couplings (d=7)
MU2=0.113; L=2
OUT="${OUT:-scan_2t_grid}"; mkdir -p "$OUT"   # GH_GUPD (if exported) propagates to gh_string -> cured run
CSV="$OUT/2t_grid.csv"
echo "beta,kappa,start,plaq,acc,Llink" > "$CSV"
BETAS="0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0"
KAPPAS="0 0.5 1.0 1.5 2.0 3.0 4.0 6.0 8.0"
NTHERM="${NTHERM:-300}"; NMEAS="${NMEAS:-100}"; NMD="${NMD:-40}"; SEED="${SEED:-12345}"

run_one() {
  b="$1"; k="$2"; st="$3"
  E="GH_FROZEN=1"; [ "$st" = cold ] && E="$E GH_COLD=1"
  o=$(env $E OMP_NUM_THREADS=2 "$BIN" 6 "$L" "$b" "$k" "$MU2" "$C" "$NTHERM" "$NMEAS" "$NMD" 1.0 "$SEED" 2>/dev/null)
  p=$(echo "$o" | grep -m1 'avg_plaquette =' | grep -oE '[-0-9.]+' | head -1)
  a=$(echo "$o" | grep -m1 'acceptance'      | grep -oE '[0-9.]+'  | head -1)
  ll=$(echo "$o"| grep -m1 'L_link'          | grep -oE '[-0-9.]+' | head -1)
  echo "$b,$k,$st,${p:-NA},${a:-NA},${ll:-NA}" >> "$CSV"
}
export -f run_one; export BIN C MU2 L OUT CSV NTHERM NMEAS NMD SEED
echo "START $(date)" > "$OUT/scan.log"
for b in $BETAS; do for k in $KAPPAS; do for st in hot cold; do echo "$b $k $st"; done; done; done \
  | xargs -P 24 -n 3 bash -c 'run_one "$@"' _
echo "ALL DONE $(date)" >> "$OUT/scan.log"
