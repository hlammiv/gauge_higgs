#!/usr/bin/env bash
# u1_scan_campaign.sh -- orchestrate a charge-q U(1) (beta,kappa) phase-diagram /
# triple-point campaign on a single workstation by fanning the grid out as many
# INDEPENDENT single-point `u1_scan` processes (the embarrassingly-parallel axis),
# K at a time, each OpenMP-threaded.
#
# WHY single-point jobs (one per L,q,beta,kappa):
#   - u1_scan seeds each node as Rng::key(base_seed, node) with node = index WITHIN
#     that invocation, and opens summary.csv with "w". Launching per-point with a
#     distinct, deterministic base_seed and its OWN output dir therefore gives
#     (a) slice-independent reproducibility, (b) no summary.csv clobbering,
#     (c) trivial resume (skip points whose output already exists), and
#     (d) the best load-balancing / core utilisation across the grid.
#
# SAFETY: DRY-RUN by default -- prints the plan, point count, and total trajectory
# budget, and launches NOTHING. Pass --run to actually execute. (These are the
# "bigger runs" that are deferred until explicitly green-lit.)
#
# Usage:
#   scripts/u1_scan_campaign.sh [--config FILE] [--out DIR] [--jobs N] [--threads N]
#                               [--run] [--force] [--build] [--aggregate] [-h]
#
# Configure the grid by editing the CONFIG block below, or pass --config FILE that
# sets any of: LS, QS, BETAS, KAPPAS (bash arrays), LAMBDA, NTHERM, NMEAS, NMD, TAU,
# MEASURE_EVERY, CAMPAIGN_SEED, NDIM, AUTOTUNE, N_SCALAR, PRETUNE, PRETUNE_SAFETY. Example config:
#     LS=(8 12); QS=(1 2 3 4 5 6 8)
#     BETAS=($(seq 0.6 0.05 1.4)); KAPPAS=($(seq 0.1 0.05 0.6))
#     LAMBDA=0.5; NTHERM=2000; NMEAS=4000; NMD=6; MEASURE_EVERY=2; CAMPAIGN_SEED=20260529
#     N_SCALAR=6; PRETUNE=1
set -euo pipefail

# ----------------------------- CONFIG (defaults) -----------------------------
# Small example grid by design, so an accidental --run does not launch a monster.
# Override via --config FILE or by editing here. q<=4 -> finite triple point;
# q>=5 -> Coulomb wedge reaches kappa=inf (no finite triple point) -- the headline.
LS=(8)                         # lattice sizes L (L^NDIM)
QS=(1 2 3 4 5)                 # Higgs charges q (= N)
BETAS=(0.8 0.9 1.0 1.1 1.2)    # gauge couplings
KAPPAS=(0.2 0.3 0.4 0.5)       # hopping couplings
LAMBDA=0.5                     # quartic (fixed across the campaign; report it)
NTHERM=1000
NMEAS=2000
NMD=20                         # nmd (start value; final value set by autotune/pretune)
TAU=1.0
MEASURE_EVERY=2
CAMPAIGN_SEED=20260529         # campaign-wide seed offset (reproducibility)
NDIM=4                         # the triple point is a D=4 object
AUTOTUNE=1                     # per-point nmd auto-tune (1=on, 0=off); IGNORED when PRETUNE=1
N_SCALAR=1                     # multi-timescale scalar sub-steps/gauge step (1=single, >1=Sexton-Weingarten)
PRETUNE=0                      # 1: pre-tune nmd ONCE per (L,q) at the stiff grid corner and reuse it
                               #    across that (beta,kappa) grid with AUTOTUNE forced 0 (recommended for
                               #    multi-timescale runs -- avoids the per-point calibration-cost blowup).
PRETUNE_NTHERM=150             # thermalization trajectories before each pre-tune burst
PRETUNE_SAFETY=1.15            # inflate the max-over-probes nmd by this factor (ceil) for headroom

# ----------------------------- runtime defaults ------------------------------
NCORES="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
THREADS=2                       # OpenMP threads per job
JOBS=""                         # concurrent jobs (default: NCORES/THREADS)
OUTROOT="u1scan_campaign"
RUN=0; FORCE=0; BUILD=0; AGG_ONLY=0
CONFIG=""

# ------------------------------- arg parsing ---------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)    CONFIG="$2"; shift 2;;
    --out)       OUTROOT="$2"; shift 2;;
    --jobs)      JOBS="$2"; shift 2;;
    --threads)   THREADS="$2"; shift 2;;
    --run)       RUN=1; shift;;
    --force)     FORCE=1; shift;;
    --build)     BUILD=1; shift;;
    --aggregate) AGG_ONLY=1; shift;;
    -h|--help)   sed -n '2,40p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[[ -n "$CONFIG" ]] && { echo "# sourcing config: $CONFIG"; source "$CONFIG"; }
[[ -z "$JOBS" ]] && JOBS=$(( NCORES / THREADS )); (( JOBS < 1 )) && JOBS=1

# repo root = parent of this script's dir
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="$REPO/build/u1_scan"

# ------------------------------ aggregation ----------------------------------
aggregate() {
  local out="$OUTROOT/summary_all.csv"; local wrote_header=0
  : > "$out"
  while IFS= read -r -d '' f; do
    if (( wrote_header == 0 )); then grep -E '^beta,' "$f" >> "$out" && wrote_header=1; fi
    grep -E '^[0-9]' "$f" >> "$out" || true
  done < <(find "$OUTROOT" -name summary.csv -not -path "$out" -print0 | sort -z)
  local n; n=$(grep -cE '^[0-9]' "$out" 2>/dev/null || echo 0)
  echo "# aggregated $n point-rows -> $out"
}
if (( AGG_ONLY )); then aggregate; exit 0; fi

# --------------------------------- build -------------------------------------
if (( BUILD )) || [[ ! -x "$BIN" ]]; then
  echo "# building u1_scan (NDIM=$NDIM) ..."
  ( cd "$REPO" && make NDIM="$NDIM" NCOL=2 build/u1_scan )
fi
[[ -x "$BIN" ]] || { echo "ERROR: $BIN not built" >&2; exit 1; }

# --------------------------- enumerate the work list -------------------------
# Canonical order L -> q -> beta -> kappa fixes the per-point global index (hence
# the seed), so results are reproducible and independent of how the grid is sliced.
declare -a JOB_L JOB_Q JOB_B JOB_K JOB_SEED JOB_DIR
idx=0; n_total=0; n_todo=0
for L in "${LS[@]}"; do
 for q in "${QS[@]}"; do
  for b in "${BETAS[@]}"; do
   for k in "${KAPPAS[@]}"; do
     bdir="$(printf 'b%.6f_k%.6f' "$b" "$k")"
     pdir="$OUTROOT/L${L}/q${q}/${bdir}"
     # per-point base_seed: campaign block in high bits, global index in low bits
     seed=$(( (CAMPAIGN_SEED << 20) + idx ))
     n_total=$((n_total+1)); idx=$((idx+1))
     # resume: a completed point has summary.csv with a numeric data row
     if (( FORCE == 0 )) && [[ -f "$pdir/summary.csv" ]] && grep -qE '^[0-9]' "$pdir/summary.csv" 2>/dev/null; then
       continue
     fi
     JOB_L+=("$L"); JOB_Q+=("$q"); JOB_B+=("$b"); JOB_K+=("$k"); JOB_SEED+=("$seed"); JOB_DIR+=("$pdir")
     n_todo=$((n_todo+1))
   done
  done
 done
done

traj_per_pt=$(( NTHERM + NMEAS * MEASURE_EVERY ))
echo "==================== U(1) scan campaign ===================="
echo "  grid: L={${LS[*]}}  q={${QS[*]}}  #beta=${#BETAS[@]}  #kappa=${#KAPPAS[@]}  -> ${n_total} points"
echo "  per point: lambda=$LAMBDA ntherm=$NTHERM nmeas=$NMEAS measure_every=$MEASURE_EVERY tau=$TAU n_scalar=$N_SCALAR  (~${traj_per_pt} traj)"
if (( PRETUNE )); then
  echo "  nmd: PRE-TUNED once per (L,q) at the stiff grid corner, reused across the grid (autotune OFF)"
else
  echo "  nmd=$NMD  autotune=$AUTOTUNE"
fi
echo "  to run: ${n_todo}  (skipped ${n_todo}!=... -> $((n_total - n_todo)) already complete)"
echo "  total trajectory budget (to-run): $(( n_todo * traj_per_pt ))"
echo "  parallelism: JOBS=$JOBS x THREADS=$THREADS  (NCORES=$NCORES)  outroot=$OUTROOT"
(( JOBS * THREADS > NCORES )) && echo "  WARNING: JOBS*THREADS=$((JOBS*THREADS)) > NCORES=$NCORES (oversubscribed)"
echo "============================================================"

if (( RUN == 0 )); then
  echo "DRY-RUN (no jobs launched). Re-run with --run to execute. Sample commands:"
  nmd_disp="$NMD"; at_disp="$AUTOTUNE"
  (( PRETUNE )) && { nmd_disp="<pretuned/Lq>"; at_disp=0; }
  for i in 0 1 2; do
    (( i >= n_todo )) && break
    echo "  OMP_NUM_THREADS=$THREADS $BIN ${JOB_L[i]} ${JOB_B[i]} ${JOB_B[i]} 1 ${JOB_K[i]} ${JOB_K[i]} 1 $LAMBDA ${JOB_Q[i]} $NTHERM $NMEAS $nmd_disp $TAU ${JOB_SEED[i]} $MEASURE_EVERY ${JOB_DIR[i]} $at_disp $N_SCALAR"
  done
  (( PRETUNE )) && echo "  (+ one u1_pretune burst per (L,q) at the stiff corner before production)"
  exit 0
fi

# ------------------------------- manifest ------------------------------------
mkdir -p "$OUTROOT"
man="$OUTROOT/manifest.txt"
{
  echo "# u1_scan campaign manifest  (NDIM=$NDIM)"
  echo "# CAMPAIGN_SEED=$CAMPAIGN_SEED LAMBDA=$LAMBDA NTHERM=$NTHERM NMEAS=$NMEAS NMD=$NMD TAU=$TAU MEASURE_EVERY=$MEASURE_EVERY AUTOTUNE=$AUTOTUNE N_SCALAR=$N_SCALAR PRETUNE=$PRETUNE"
  echo "# columns: L q beta kappa base_seed outdir"
  for i in "${!JOB_L[@]}"; do
    echo "${JOB_L[i]} ${JOB_Q[i]} ${JOB_B[i]} ${JOB_K[i]} ${JOB_SEED[i]} ${JOB_DIR[i]}"
  done
} > "$man"
echo "# manifest -> $man"

# ----------------------- per-(L,q) nmd PRE-TUNE ------------------------------
# Tune nmd ONCE per (L,q) and reuse it across the grid with autotune OFF (avoids
# the per-point calibration-cost blowup). The multi-timescale integrator makes nmd
# ~flat in (beta,kappa), but NOT perfectly: the MD forces are stiffest near the
# TRANSITION (interior), and a deep-Higgs corner at max-kappa can be SMOOTHER than
# the transition row just below it. A corner-only tune therefore UNDERESTIMATES nmd
# and can starve a transition row (the q5 kappa=0.9 acc~3e-4 blowup). So we probe a
# SPREAD (max-kappa corners + interior midpoints), take the MAX nmd, and inflate by
# PRETUNE_SAFETY for headroom. Cost: |probes| pretune bursts per (L,q) (~|probes| x
# (PRETUNE_NTHERM+30) traj) -- cheap vs production, but scales with L^NDIM.
declare -A NMD_LQ
if (( PRETUNE )); then
  PBIN="$REPO/build/u1_pretune"
  if [[ ! -x "$PBIN" ]]; then
    echo "# building u1_pretune (NDIM=$NDIM) ..."
    ( cd "$REPO" && make NDIM="$NDIM" NCOL=2 build/u1_pretune )
  fi
  [[ -x "$PBIN" ]] || { echo "ERROR: $PBIN not built" >&2; exit 1; }
  mapfile -t _bs < <(printf '%s\n' "${BETAS[@]}"  | sort -g)
  mapfile -t _ks < <(printf '%s\n' "${KAPPAS[@]}" | sort -g)
  bmin="${_bs[0]}"; bmax="${_bs[$(( ${#_bs[@]} - 1 ))]}"; bmid="${_bs[$(( ${#_bs[@]} / 2 ))]}"
  kmid="${_ks[$(( ${#_ks[@]} / 2 ))]}"; kmax="${_ks[$(( ${#_ks[@]} - 1 ))]}"
  # Probe set: both max-kappa corners + interior of the top row + mid-plane spread.
  PROBES=( "$bmax $kmax" "$bmin $kmax" "$bmid $kmax" "$bmid $kmid" "$bmin $kmid" "$bmax $kmid" )
  plog="$OUTROOT/pretune.log"; : > "$plog"
  echo "# pre-tuning nmd per (L,q) over ${#PROBES[@]} probes (corners + interior), max x safety=$PRETUNE_SAFETY (diagnostics -> $plog) ..."
  { echo "# pretune probes (beta kappa): ${PROBES[*]}  | lambda=$LAMBDA n_scalar=$N_SCALAR ntherm=$PRETUNE_NTHERM safety=$PRETUNE_SAFETY"; } >> "$plog"
  for L in "${LS[@]}"; do for q in "${QS[@]}"; do
    pbase=$(( (CAMPAIGN_SEED << 20) + 900000 + L*100 + q ))
    nmd=0; detail=""; pj=0
    for probe in "${PROBES[@]}"; do
      read -r pb pk <<< "$probe"
      ni=$(OMP_NUM_THREADS="$THREADS" "$PBIN" "$L" "$pb" "$pk" "$LAMBDA" "$q" "$PRETUNE_NTHERM" "$NMD" "$TAU" "$(( pbase + pj ))" "$N_SCALAR" 2>>"$plog")
      [[ "$ni" =~ ^[0-9]+$ ]] || ni="$NMD"
      (( ni > nmd )) && nmd="$ni"
      detail+=" (b$pb,k$pk)=$ni"
      pj=$((pj+1))
    done
    # safety inflation (ceil); fall back to NMD only if the probes all failed.
    nmd=$(awk -v n="$nmd" -v s="$PRETUNE_SAFETY" 'BEGIN{ v=n*s; r=int(v); if (v>r) r++; print r }')
    if (( nmd < 1 )); then nmd="$NMD"; fi
    NMD_LQ["${L}_${q}"]="$nmd"
    echo "#   (L=$L, q=$q) -> nmd=$nmd   [max-over-probes x $PRETUNE_SAFETY;$detail]" | tee -a "$plog"
  done; done
fi

# --------------------------------- launch ------------------------------------
run_point() {  # L b k q seed dir
  local L="$1" b="$2" k="$3" q="$4" seed="$5" dir="$6"
  local nmd_use="$NMD" at_use="$AUTOTUNE"
  if (( PRETUNE )); then nmd_use="${NMD_LQ[${L}_${q}]:-}"; at_use=0; fi
  if ! [[ "$nmd_use" =~ ^[0-9]+$ ]]; then nmd_use="$NMD"; fi
  mkdir -p "$dir"
  OMP_NUM_THREADS="$THREADS" "$BIN" "$L" "$b" "$b" 1 "$k" "$k" 1 "$LAMBDA" "$q" \
      "$NTHERM" "$NMEAS" "$nmd_use" "$TAU" "$seed" "$MEASURE_EVERY" "$dir" "$at_use" "$N_SCALAR" \
      > "$dir/run.log" 2>&1
}

echo "# launching $n_todo jobs, $JOBS at a time ..."
done_ct=0
for i in "${!JOB_L[@]}"; do
  while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n 2>/dev/null || sleep 1; done
  run_point "${JOB_L[i]}" "${JOB_B[i]}" "${JOB_K[i]}" "${JOB_Q[i]}" "${JOB_SEED[i]}" "${JOB_DIR[i]}" &
done
wait
echo "# all jobs finished."
aggregate
echo "# DONE. Per-point time series + summary.csv under $OUTROOT/ ; combined -> $OUTROOT/summary_all.csv"
echo "# next: feed the per-point ts_*.dat (columns: traj A B ...) to reweight.hpp / autocorr.hpp."
