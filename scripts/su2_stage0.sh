#!/usr/bin/env bash
# SU(2)->H Stage-0 shakedown: HMC health (acceptance, <exp(-dH)>~1), per-traj COST, and
# the L_link low-vs-high-kappa CONTRAST (condensation responds) for the 5 reps:
#   adjoint->U(1) (control) / spin-2->Q8 soft (control) / 2T / 2O / 2I.
# At 2^4 and 4^4 (2I: 2^4 only -- 56 s/traj at 4^4 at OMP=1). All OMP=1, run concurrently
# (18 points <= cores). NOT physics-grade stats -- a shakedown to size the real scan.
set -uo pipefail
cd "$(dirname "$0")/.."
BIN=build/hmc_higgs_multi
[[ -x $BIN ]] || { echo "ERROR: build $BIN first (make NDIM=4 NCOL=2 build/hmc_higgs_multi)"; exit 1; }
OUT=su2_stage0; mkdir -p "$OUT"
SUM="$OUT/summary.txt"; : > "$SUM"
BETA=2.0; NTHERM=80; NMEAS=120; NMD=20; TAU=1.0
KAPPAS=(0.1 0.6)

# rep : label : mu2 : couplings(comma list, or 'auto')   [f_c from docs/locking_couplings.md]
SPECS=(
  "adj:U1adj:0.1:auto"
  "4:Q8soft:0.2225:0.2964,0.2648,0.0970,0.2946,0.0472"
  "6:2T:0.113:0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
  "8:2O:0.108:0.0562,0.2045,0.2230,0.2206,0.0444,0.0740,0.0066,0.0340,0.1366"
  "12:2I:0.065:0.1182,0.0819,0.0353,0.0274,0.1550,0.0828,0.0157,0.0702,0.1083,0.1143,0.0189,0.1124,0.0595"
)

run_one() {
  local rep=$1 lab=$2 mu2=$3 coup=$4 L=$5 kap=$6 seed=$7
  local log="$OUT/${lab}_L${L}_k${kap}.log"
  local t0 t1 secs
  t0=$(date +%s.%N)
  OMP_NUM_THREADS=1 $BIN "$rep" "$L" "$BETA" "$kap" "$mu2" "$coup" "$NTHERM" "$NMEAS" "$NMD" "$TAU" "$seed" > "$log" 2>&1
  t1=$(date +%s.%N)
  secs=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f", b-a}')
  local plaq lphi llink acc edh
  plaq=$(awk '/^plaquette/{print $3}' "$log")
  lphi=$(awk '/^L_phi/{print $3}' "$log")
  llink=$(awk '/^L_link/{print $3}' "$log")
  acc=$(awk '/^acceptance/{print $3}' "$log")
  edh=$(awk -F'=' '/^acceptance/{print $NF}' "$log")
  printf "%-7s L=%s k=%-3s | plaq=%-8s Lphi=%-9s Llink=%-9s acc=%-6s <e-dH>=%-8s | %ss/%dtraj\n" \
    "$lab" "$L" "$kap" "${plaq:-NA}" "${lphi:-NA}" "${llink:-NA}" "${acc:-NA}" "${edh:-NA}" "$secs" "$((NTHERM+NMEAS))" >> "$SUM"
}

i=0
for spec in "${SPECS[@]}"; do
  IFS=':' read -r rep lab mu2 coup <<< "$spec"
  for L in 2 4; do
    [[ "$rep" == "12" && "$L" == "4" ]] && continue   # skip 2I 4^4 (56 s/traj at OMP=1)
    for kap in "${KAPPAS[@]}"; do
      i=$((i+1)); run_one "$rep" "$lab" "$mu2" "$coup" "$L" "$kap" "$((1000+i))" &
    done
  done
done
wait
echo "==================== SU(2) Stage-0 shakedown (beta=$BETA, ntherm=$NTHERM nmeas=$NMEAS nmd=$NMD) ===================="
sort "$SUM"
echo "================================================================================"
echo "health: want acc in ~0.6-0.95 and <e-dH> ~ 1.0 ; contrast: Llink(k=0.6) > Llink(k=0.1)"
