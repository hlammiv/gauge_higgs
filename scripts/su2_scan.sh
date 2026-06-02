#!/usr/bin/env bash
# SU(2)->H WIDE (beta,kappa) phase scan at 2^4 -- decorrelation-instrumented + distributable.
# Uniform 2^4 so all 5 reps (adjoint U(1) / Q8 soft / 2T / 2O / 2I) are volume-comparable.
# Each point writes a per-trajectory time series (OUT/ts/) and a CSV row INCLUDING the
# driver's self-reported tau_int(L_link) + N_eff, so decorrelation is verified everywhere.
# Adaptive nmd (rises with kappa; 2I gets a stiff explicit table). The instrumented
# build/hmc_higgs_multi reports tau_int and uses the autocorrelation-aware L_link error.
#
# DISTRIBUTABLE via env vars (same script on local / lenore / loranne, disjoint beta-slices):
#   REPS="2I"  BETAS="2.5 2.75 3.0"  JOBS=4  OUT=su2scan_w  bash scripts/su2_scan.sh
# Env-overridable: REPS, BETAS, KAPPAS, NTHERM, NMEAS, JOBS, OUT, SEED0, BIN.
set -uo pipefail
cd "$(dirname "$0")/.."
BIN=${BIN:-build/hmc_higgs_multi}
[[ -x $BIN ]] || { echo "ERROR: $BIN not built (make NDIM=4 NCOL=2 build/hmc_higgs_multi)"; exit 1; }
OUT=${OUT:-su2scan_w}; mkdir -p "$OUT" "$OUT/ts"
NTHERM=${NTHERM:-500}; NMEAS=${NMEAS:-800}; TAU=1.0; SEED0=${SEED0:-20260602}; L=${L:-2}
# GH_FROZEN (env) is inherited by the driver subprocesses -> |phi|=1 frozen-length scalar.
# Unquoted assignment word-splits on IFS (space AND newline), so this handles BOTH a
# space-separated env override (BETAS="2.5 2.75 3.0") and seq's newline-separated default.
# shellcheck disable=SC2206
BETAS=(  ${BETAS:-$(seq -f '%.2f' 1.0 0.25 3.0)}  )     # 9: 1.00..3.00
KAPPAS=( ${KAPPAS:-$(seq -f '%.1f' 0.0 0.1 1.0)}  )     # 11: 0.0..1.0
REPS=(   ${REPS:-adj Q8soft 2T 2O 2I}             )
JOBS=${JOBS:-$(nproc)}
# Optional kappa cap (space-free, survives nested ssh): e.g. KMAX=0.7 for the costly 2I.
if [[ -n "${KMAX:-}" ]]; then
  _kf=(); for k in "${KAPPAS[@]}"; do awk -v a="$k" -v b="$KMAX" 'BEGIN{exit !(a<=b+1e-9)}' && _kf+=("$k"); done
  KAPPAS=("${_kf[@]}")
fi

# label -> "driverrep:mu2:couplings:basenmd"  (f_c from docs/locking_couplings.md)
declare -A SPEC
SPEC[adj]="adj:0.1:auto:16"
SPEC[Q8soft]="4:0.2225:0.2964,0.2648,0.0970,0.2946,0.0472:20"
SPEC[2T]="6:0.113:0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130:28"
SPEC[2O]="8:0.108:0.0562,0.2045,0.2230,0.2206,0.0444,0.0740,0.0066,0.0340,0.1366:36"
SPEC[2I]="12:0.065:0.1182,0.0819,0.0353,0.0274,0.1550,0.0828,0.0157,0.0702,0.1083,0.1143,0.0189,0.1124,0.0595:0"

# adaptive nmd(label,kappa,base): 2I uses an explicit stiff table (d=13 needs fine steps deep
# in the Higgs phase); other reps scale base up by 60% above kappa~0.5.
nmd_for() {
  if [[ "$1" == "2I" ]]; then
    awk -v k="$2" 'BEGIN{ if(k<=0.21)print 16; else if(k<=0.31)print 28; else if(k<=0.51)print 44; else if(k<=0.71)print 56; else print 68 }'
  else
    awk -v k="$2" -v b="$3" 'BEGIN{ f=(k>0.51)?1.6:1.0; print int(b*f+0.999) }'
  fi
}

run_one() {  # label drrep mu2 coup base beta kappa seed
  local lab=$1 dr=$2 mu2=$3 coup=$4 base=$5 beta=$6 kap=$7 seed=$8
  local nmd; nmd=$(nmd_for "$lab" "$kap" "$base")
  local log="$OUT/${lab}_b${beta}_k${kap}.log" ts="$OUT/ts/${lab}_b${beta}_k${kap}.dat"
  OMP_NUM_THREADS=1 $BIN "$dr" "$L" "$beta" "$kap" "$mu2" "$coup" "$NTHERM" "$NMEAS" "$nmd" "$TAU" "$seed" "$ts" > "$log" 2>&1
  awk -v b="$beta" -v k="$kap" -v nm="$nmd" '
    /^plaquette/{p=$3;pe=$5}
    /^L_phi/    {lp=$3}
    /^L_link/   {ll=$3;le=$5;for(i=1;i<=NF;i++)if($i~/^chi_link=/){split($i,a,"=");ch=a[2]}}
    /^autocorr/ {for(i=1;i<=NF;i++){if($i~/tau_int/){split($i,t,"=");ti=t[2]} if($i~/^N_eff=/){split($i,n,"=");ne=n[2]}}}
    /^acceptance/{ac=$3;m=split($0,q,"=");edh=q[m]}
    END{printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",b,k,nm,p,pe,lp,ll,le,ch,ti,ne,ac,edh}' "$log" >> "$OUT/${lab}.csv"
}

# headers (only create if absent, so distributed slices append to the same CSV safely)
for lab in "${REPS[@]}"; do
  f="$OUT/${lab}.csv"
  [[ -f $f ]] || echo "beta,kappa,nmd,plaq,plaq_err,Lphi,Llink,Llink_err,chi_link,tau_int,N_eff,acc,edh" > "$f"
done

echo "# WIDE scan: reps={${REPS[*]}} beta={${BETAS[*]}} kappa={${KAPPAS[*]}} -> $(( ${#REPS[@]} * ${#BETAS[@]} * ${#KAPPAS[@]} )) pts"
echo "# ntherm=$NTHERM nmeas=$NMEAS L=$L^4 JOBS=$JOBS OUT=$OUT (adaptive nmd; ts in $OUT/ts/)"
i=0
for lab in "${REPS[@]}"; do
  IFS=':' read -r dr mu2 coup base <<< "${SPEC[$lab]}"
  for beta in "${BETAS[@]}"; do for kap in "${KAPPAS[@]}"; do
    while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n 2>/dev/null || sleep 1; done
    i=$((i+1)); run_one "$lab" "$dr" "$mu2" "$coup" "$base" "$beta" "$kap" "$((SEED0+i))" &
  done; done
done
wait
echo "# DONE reps={${REPS[*]}}:"
for lab in "${REPS[@]}"; do echo "  $OUT/${lab}.csv ($(( $(wc -l < "$OUT/${lab}.csv") - 1 )) rows)"; done
