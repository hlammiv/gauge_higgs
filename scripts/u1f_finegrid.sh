#!/usr/bin/env bash
# FINE-GRID campaign for lenore + lucia. Two goals:
#   (1) the TWO Coulomb metrics agree everywhere -> combined `pm` run (Rmax>1) emits m_gamma AND the spatial
#       Wilson V(R) from the SAME config stream, so both come from identical statistics/lattice (no geometry
#       confound). Every deconfined point then carries both metrics.
#   (2) denser grid near the transitions -> two tiers:
#       TIER A (WALL): cheap `point` (L=12) on a FINE (beta,kappa) mesh through the confinement wall +
#                      triple point -> sharp sigma_1 / <cos> / rho_M lines.
#       TIER B (COUL): combined `pm`+V(R) (L_s=20, Lt=8, Rmax=10) on the deconfined/wedge region -> both
#                      Coulomb metrics, denser in beta near the wall.
# Jobs are split across the two machines (core-weighted), WALL first (fast) then COUL. Resumable (# DONE skip).
#
# *** DRY RUN BY DEFAULT *** : prints the plan + the exact launch commands and writes the queue files, but
# launches NOTHING. Re-run with GO=1 to actually sync/build/launch.
#   inspect:  scripts/u1f_finegrid.sh
#   launch :  GO=1 scripts/u1f_finegrid.sh
# tune via env (examples): QS="4 5 6 8"  COUL_KAPS="1.0 1.5 2.5"  WALL_DB=0.1   OUT=u1f_fg
set -u
OUT="${OUT:-u1f_fg}"
LEN_CORES="${LEN_CORES:-30}"; LUC_CORES="${LUC_CORES:-18}"
QS="${QS:-2 3 4 5 6 8}"
# Tier A (wall/triple-point): point mode, cheap, dense
WALL_BETAS="${WALL_BETAS:-0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45}"
WALL_KAPS="${WALL_KAPS:-0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9}"
WALL_L="${WALL_L:-12}"; WALL_NS="${WALL_NS:-4000}"; WALL_NT="${WALL_NT:-1500}"
# Tier B (Coulomb region): combined pm+V(R), denser in beta near the wall
COUL_BETAS="${COUL_BETAS:-1.1 1.2 1.3 1.4 1.5 1.7 2.0 2.5}"
COUL_KAPS="${COUL_KAPS:-0.6 1.0 1.5 2.0 2.5}"
COUL_LS="${COUL_LS:-20}"; COUL_LT="${COUL_LT:-8}"; COUL_RMAX="${COUL_RMAX:-8}"
COUL_NS="${COUL_NS:-3000}"; COUL_NT="${COUL_NT:-1500}"; COUL_ME="${COUL_ME:-5}"
mkdir -p "$OUT"
QL="$OUT/queue_lenore.q"; QU="$OUT/queue_lucia.q"; : > "$QL"; : > "$QU"

# core-weighted round-robin: of every (LEN+LUC) jobs, the first LEN_CORES go to lenore, the rest to lucia
emit() {   # $1 = job line
  local tot=$((LEN_CORES + LUC_CORES))
  if [ $(( I % tot )) -lt $LEN_CORES ]; then echo "$1" >> "$QL"; else echo "$1" >> "$QU"; fi
  I=$((I + 1))
}

I=0; nwall=0
for q in $QS; do for b in $WALL_BETAS; do for k in $WALL_KAPS; do
  seed=$((300000 + I))
  emit "point $WALL_L $b $k $q $WALL_NS $WALL_NT $seed 3"; nwall=$((nwall + 1))
done; done; done
ncoul=0
for q in $QS; do for b in $COUL_BETAS; do for k in $COUL_KAPS; do
  seed=$((400000 + I))
  emit "pm $COUL_LS $COUL_LT $b $k $q $COUL_NS $COUL_NT $COUL_ME $seed 3 $COUL_RMAX"; ncoul=$((ncoul + 1))
done; done; done

NL=$(wc -l < "$QL"); NU=$(wc -l < "$QU")
cat <<EOF
=== fine-grid plan (OUT=$OUT) ===
TIER A (wall/triple point): point L=$WALL_L  betas[$(echo $WALL_BETAS|wc -w)] x kaps[$(echo $WALL_KAPS|wc -w)] x q[$(echo $QS|wc -w)] = $nwall jobs  (~1-2 min each)
TIER B (Coulomb metrics) : pm+V(R) L_s=$COUL_LS Lt=$COUL_LT Rmax=$COUL_RMAX  betas[$(echo $COUL_BETAS|wc -w)] x kaps[$(echo $COUL_KAPS|wc -w)] x q[$(echo $QS|wc -w)] = $ncoul jobs  (~60-90 min each, single-thread)
split: lenore=$NL jobs ($LEN_CORES cores) | lucia=$NU jobs ($LUC_CORES cores)
rough wall-clock: TIER A ~$(( (nwall/(LEN_CORES+LUC_CORES)+1)*2 )) min ; TIER B ~$(( (ncoul/(LEN_CORES+LUC_CORES)+1)*75 )) min
queues written: $QL  $QU
EOF

if [ "${GO:-0}" != "1" ]; then
  echo
  echo "DRY RUN -- nothing launched. To launch:  GO=1 $0"
  echo "would run, lucia (local):   POT_OUT=$OUT scripts/u1f_launch.sh $QU $LUC_CORES $OUT"
  echo "would run, lenore (remote): rsync src+queue -> build -> scripts/u1f_launch.sh queue_lenore $LEN_CORES $OUT"
  exit 0
fi

echo "=== LAUNCHING ==="
# lucia (local)
scripts/u1f_launch.sh "$QU" "$LUC_CORES" "$OUT"
# lenore (remote): sync source + queue, rebuild, launch
rsync -az -e 'ssh -p 60022' src/ lenore_remote:higgs_gauge/src/ >/dev/null 2>&1
rsync -az -e 'ssh -p 60022' scripts/u1f_job.sh scripts/u1f_launch.sh "$QL" lenore_remote:higgs_gauge/scripts/ >/dev/null 2>&1
ssh -p 60022 lenore_remote "cd ~/higgs_gauge && g++ -O3 -march=native -std=c++17 -fopenmp -Isrc src/u1_frozen.cpp -o build/u1_frozen && scripts/u1f_launch.sh scripts/$(basename "$QL") $LEN_CORES $OUT"
echo "launched on both. progress: grep -lc '# DONE' $OUT/job_*.out | wc -l  of  $((NL+NU))"
