#!/usr/bin/env bash
# Benchmark the GaugeHiggsHMC integrator: cost-per-accepted-trajectory vs nmd, at a
# chosen (rep, L, beta, kappa). Sweeps nmd to map the acceptance/cost frontier so an
# integrator change (e.g. force-gradient) can be compared HEAD-TO-HEAD at matched
# acceptance, not at fixed nmd (which is misleading).
#
# Usage: gh_integrator_bench.sh <binary> <rep> <L> <beta> <kappa> <couplings> "<nmd_list>" [extra_env]
#   binary    : path to hmc_higgs_multi (production or a scratch prototype build)
#   extra_env  : e.g. "GH_FG=1" to select a prototype integrator branch
# Always runs frozen (GH_FROZEN=1), OMP_NUM_THREADS=1, fixed seed=7.
#
# The winning integrator is the one with the LOWEST cost/acc (= s_per_traj / acceptance)
# at the target acceptance band (~0.65-0.85). Reaching higher acceptance at higher nmd is
# NOT a win if the cost/acc is worse. Force-eval count = nmd * (2 for 2MN, 3 for FG).
set -u
BIN=${1:?binary}; REP=${2:?rep}; L=${3:?L}; BETA=${4:?beta}; KAPPA=${5:?kappa}
CC=${6:?couplings}; NMDS=${7:?nmd list}; EENV=${8:-}
NTH=${NTH:-15}; NME=${NME:-60}
printf "# bench bin=%s rep=%s L=%s beta=%s kappa=%s env='%s'  nth=%s nme=%s\n" \
  "$BIN" "$REP" "$L" "$BETA" "$KAPPA" "$EENV" "$NTH" "$NME"
printf "# %-6s %-12s %-10s %-10s\n" nmd s/traj acc cost/acc
for nmd in $NMDS; do
  t0=$(date +%s.%N)
  out=$(eval "$EENV GH_FROZEN=1 OMP_NUM_THREADS=1 $BIN $REP $L $BETA $KAPPA 0.113 \"$CC\" $NTH $NME $nmd 1 7" 2>&1)
  t1=$(date +%s.%N)
  acc=$(echo "$out" | grep -oP 'acceptance\s*=\s*\K[0-9.]+')
  dt=$(echo "$t1 - $t0" | bc); ntot=$((NTH+NME)); per=$(echo "scale=4; $dt/$ntot" | bc)
  cpa=$(echo "scale=4; $per/(${acc:-0}+0.0001)" | bc)
  printf "  %-6s %-12s %-10s %-10s\n" "$nmd" "$per" "${acc:-NA}" "$cpa"
done
