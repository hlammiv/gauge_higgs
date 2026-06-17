#!/usr/bin/env bash
# Benchmark + bit-identity harness for the "share the forward D(U)phi_y between the scalar-
# hopping force and the matter-link back-reaction" optimization (GaugeHiggsHMC tensor path,
# BT spin-3 {6} / BO spin-4 {8}).  Builds a BASELINE driver from the given production src and
# a DEDUP driver from a prototype src, then (1) checks bit-identity of observables at matched
# seed/params and (2) times s/traj across an nmd sweep, OMP=1.
#
# Usage:
#   scripts/gh_force_dedup_bench.sh <baseline_src_dir> <dedup_src_dir> [rep] [L] [ntraj]
# e.g.
#   scripts/gh_force_dedup_bench.sh src /tmp/scratch_dedup/src 6 6 40
#
# rep: 6=spin-3(d=7,BT), 8=spin-4(d=9,BO), 12=spin-6(d=13,BI -> fast path, dedup is a no-op).
set -euo pipefail
BASE_SRC=${1:?baseline src dir}; DEDUP_SRC=${2:?dedup src dir}
REP=${3:-6}; L=${4:-6}; NTRAJ=${5:-40}
CXX=${CXX:-g++}; FLAGS="-std=c++20 -O3 -march=native -funroll-loops -fopenmp -DNDIM=4 -DNCOL=2"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
echo "# building baseline ($BASE_SRC) + dedup ($DEDUP_SRC) ..."
$CXX $FLAGS -I"$BASE_SRC"  -o "$TMP/base"  "$BASE_SRC/hmc_higgs_multi.cpp"
$CXX $FLAGS -I"$DEDUP_SRC" -o "$TMP/dedup" "$DEDUP_SRC/hmc_higgs_multi.cpp"

# 7-channel couplings for spin-3; "auto" otherwise (channel count varies by rep).
COUP="0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130"
[ "$REP" = "6" ] || COUP="auto"

echo "# --- bit-identity (rep=$REP L=$L frozen, seed 12345, $NTRAJ meas) ---"
ARGS="$REP $L 1.5 1.0 0.113 $COUP 20 $NTRAJ 8 1 12345"
B=$(GH_FROZEN=1 OMP_NUM_THREADS=1 "$TMP/base"  $ARGS 2>/dev/null | grep -E 'plaquette|L_link|acceptance')
D=$(GH_FROZEN=1 OMP_NUM_THREADS=1 "$TMP/dedup" $ARGS 2>/dev/null | grep -E 'plaquette|L_link|acceptance')
if [ "$B" = "$D" ]; then echo "BIT-IDENTICAL ✓"; else echo "MISMATCH ✗"; echo "base: $B"; echo "dedup:$D"; fi

echo "# --- timing (rep=$REP L=$L frozen, OMP=1, $NTRAJ trajs) ---"
for nmd in 8 16; do
  for bin in base dedup; do
    t0=$(date +%s.%N)
    GH_FROZEN=1 OMP_NUM_THREADS=1 "$TMP/$bin" $REP $L 1.5 1.0 0.113 $COUP 0 $NTRAJ $nmd 1 1 >/dev/null 2>&1
    t1=$(date +%s.%N)
    printf "  nmd=%2d %-6s %.4f s/traj\n" "$nmd" "$bin" "$(echo "scale=5;($t1-$t0)/$NTRAJ"|bc)"
  done
done
