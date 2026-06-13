#!/usr/bin/env bash
# PHASE-3 REFINEMENT campaign (combined pm+V(R)) to sharpen the (beta,kappa) phase diagram. Two goals:
#   (1) ASYMPTOTE: extend kappa high (3.5/5/7) so the Coulomb-Higgs boundary for q=5,6 (and 4,8) flattens to
#       its kappa->inf value -- at kappa=2.5 it is still drifting.
#   (2) SMOOTH: finer beta (Delta=0.1 through the transition band) so the heatmaps are smooth at the boundary.
# Same lattice/stats as the existing u1f_fg grid (L_s=20,Lt=8,Rmax=8,nsweep=3000) so cells MERGE homogeneously.
# SKIP-EXISTING: only generates target cells not already DONE in the local u1f_fg (so no recompute). Per-cell
# deterministic seeds (900000-base) so the refinement itself resumes. Core-weighted split lenore+lucia.
#
# *** DRY RUN BY DEFAULT *** (prints plan + queues, launches nothing). GO=1 to sync/launch.
#   inspect:  scripts/u1f_phase3_refine.sh
#   priority (q5,q6 high-kappa only):  QS="5 6" KAPS="3.5 5.0 7.0" scripts/u1f_phase3_refine.sh
#   launch :  GO=1 scripts/u1f_phase3_refine.sh
set -u
OUT="${OUT:-u1f_fg}"                       # same dir as the existing grid (cells merge)
LEN_CORES="${LEN_CORES:-30}"; LUC_CORES="${LUC_CORES:-18}"
QS="${QS:-4 5 6 8}"                                                  # wedge-relevant q only (q=2,3 are flat Higgs)
BETAS="${BETAS:-1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.0 2.2 2.5}"        # finer through transition + tail
KAPS="${KAPS:-1.0 1.5 2.0 2.5 3.5 5.0 7.0}"                          # existing mid + extended high (asymptote)
LS="${LS:-20}"; LT="${LT:-8}"; RMAX="${RMAX:-8}"; NS="${NS:-3000}"; NT="${NT:-1500}"; ME="${ME:-5}"
mkdir -p "$OUT"
QL="$OUT/refine_lenore.q"; QU="$OUT/refine_lucia.q"; : > "$QL"; : > "$QU"

done_cell() {   # is (q,beta,kappa) already DONE in $OUT (any seed)?  args: q beta kappa
  local q="$1" bu="${2/./_}" ku="${3/./_}"
  local g="$OUT"/job_pm_${LS}_${LT}_${bu}_${ku}_${q}_*.out
  ls $g >/dev/null 2>&1 && grep -lq '# DONE' $g 2>/dev/null
}

I=0; new=0; skip=0
for q in $QS; do for b in $BETAS; do for k in $KAPS; do
  if done_cell "$q" "$b" "$k"; then skip=$((skip+1)); continue; fi
  seed=$((900000 + I)); I=$((I+1)); new=$((new+1))
  line="pm $LS $LT $b $k $q $NS $NT $ME $seed 3 $RMAX"
  if [ $(( (new-1) % (LEN_CORES+LUC_CORES) )) -lt $LEN_CORES ]; then echo "$line" >> "$QL"; else echo "$line" >> "$QU"; fi
done; done; done

NL=$(wc -l < "$QL"); NU=$(wc -l < "$QU"); TGT=$(( $(echo $QS|wc -w) * $(echo $BETAS|wc -w) * $(echo $KAPS|wc -w) ))
cat <<EOF
=== phase-3 refinement plan (OUT=$OUT) ===
target grid: q[$(echo $QS|wc -w)] x beta[$(echo $BETAS|wc -w)] x kappa[$(echo $KAPS|wc -w)] = $TGT cells
  beta : $BETAS
  kappa: $KAPS  (new high-kappa = those >2.5)
already DONE (skipped): $skip   |   NEW cells to run: $new
split: lenore=$NL ($LEN_CORES cores) | lucia=$NU ($LUC_CORES cores)
combined pm+V(R) ~100 min/job -> rough wall-clock ~$(( (new/(LEN_CORES+LUC_CORES)+1)*100 )) min on both machines
queues: $QL  $QU
EOF

if [ "${GO:-0}" != "1" ]; then
  echo; echo "DRY RUN -- nothing launched. Launch with:  GO=1 $0"
  echo "  (rsync lenore u1f_fg -> local FIRST so skip-existing is accurate)"
  exit 0
fi
echo "=== LAUNCHING ==="
[ "$NU" -gt 0 ] && scripts/u1f_launch.sh "$QU" "$LUC_CORES" "$OUT"
if [ "$NL" -gt 0 ]; then
  rsync -az -e 'ssh -p 60022' src/ lenore_remote:higgs_gauge/src/ >/dev/null 2>&1
  rsync -az -e 'ssh -p 60022' scripts/u1f_job.sh scripts/u1f_launch.sh "$QL" lenore_remote:higgs_gauge/scripts/ >/dev/null 2>&1
  ssh -p 60022 lenore_remote "cd ~/higgs_gauge && g++ -O3 -march=native -std=c++17 -fopenmp -Isrc src/u1_frozen.cpp -o build/u1_frozen && scripts/u1f_launch.sh scripts/$(basename "$QL") $LEN_CORES $OUT"
fi
echo "launched. progress: grep -lc '# DONE' $OUT/job_pm_*_9000*.out | wc -l  of  $new"
