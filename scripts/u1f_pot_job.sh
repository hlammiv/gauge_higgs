#!/usr/bin/env bash
# Run ONE frozen `pot` (static-potential V(R)) job single-threaded -> $POT_OUT/job_<args>.out.
# Skips if the output already completed (has the `pot` header), so a relaunch resumes the queue without
# redoing finished work. Mirrors scripts/u1f_job.sh but with the pot-mode completion marker.
set -u
OUT="${POT_OUT:-u1f_vr}"
BIN="${U1F_BIN:-./build/u1_frozen}"
mkdir -p "$OUT"
name=$(echo "$*" | tr ' /.' '___')
f="$OUT/job_${name}.out"
if [ -s "$f" ] && grep -q '# u1_frozen pot:' "$f" 2>/dev/null; then exit 0; fi
OMP_NUM_THREADS=1 "$BIN" "$@" > "$f" 2>&1
