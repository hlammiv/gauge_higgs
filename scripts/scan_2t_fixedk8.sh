#!/bin/bash
# SU(2)->2T at FIXED kappa=8 (deep Higgs / 2T-locked), beta-scan, L=2, metastability-aware:
# hot+cold starts x 2 seeds per beta. Tests whether the deep-Higgs discrete (2T) residual is
# metastable like U(1) q>=2 (hot/cold gap). Measures <plaq> + L_link (no sigma at L=2, Rmax=1).
set -u
cd ~/gh_scan
BIN=build/gh_string
C="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"   # 2T-locking quartic couplings (d=7)
MU2=0.113; KAPPA=8; L=2
OUT=scan_2t_fixedk8; mkdir -p "$OUT"
CSV="$OUT/2t_k8.csv"
echo "beta,start,seed,plaq,acc,Llink" > "$CSV"
BETAS=$(seq -f "%.2f" 0 0.25 5.0)
NTHERM=300; NMEAS=100; NMD=40

run_one() {
  b="$1"; st="$2"; seed="$3"
  E="GH_FROZEN=1"; [ "$st" = cold ] && E="$E GH_COLD=1"
  o=$(env $E OMP_NUM_THREADS=2 "$BIN" 6 "$L" "$b" "$KAPPA" "$MU2" "$C" "$NTHERM" "$NMEAS" "$NMD" 1.0 "$seed" 2>/dev/null)
  p=$(echo "$o" | grep -m1 'avg_plaquette =' | grep -oE '[-0-9.]+' | head -1)
  a=$(echo "$o" | grep -m1 'acceptance'      | grep -oE '[0-9.]+'  | head -1)
  ll=$(echo "$o"| grep -m1 'L_link'          | grep -oE '[-0-9.]+' | head -1)
  echo "$b,$st,$seed,${p:-NA},${a:-NA},${ll:-NA}" >> "$CSV"
}
export -f run_one; export BIN C MU2 KAPPA L OUT CSV NTHERM NMEAS NMD
echo "START $(date)" > "$OUT/scan.log"
for b in $BETAS; do for st in hot cold; do for seed in 12345 67890; do echo "$b $st $seed"; done; done; done \
  | xargs -P 16 -n 3 bash -c 'run_one "$@"' _
echo "ALL DONE $(date)" >> "$OUT/scan.log"
