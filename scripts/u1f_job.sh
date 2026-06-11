#!/usr/bin/env bash
# Run ONE frozen-U(1) job single-threaded; write output to $U1F_OUT/job_<args>.out.
# Args = the u1_frozen CLI (e.g. `point 8 0.9 0.3 2 4000 1500 7`). Skips if the output already
# exists and ended cleanly (so a relaunch resumes the queue without redoing finished work).
set -u
OUT="${U1F_OUT:-u1f_campaign}"
BIN="${U1F_BIN:-./build/u1_frozen}"
mkdir -p "$OUT"
name=$(echo "$*" | tr ' /.' '___')
f="$OUT/job_${name}.out"
if [ -s "$f" ] && grep -q '# DONE\|# <A>\|<cos>_down' "$f" 2>/dev/null; then exit 0; fi
OMP_NUM_THREADS=1 "$BIN" "$@" > "$f" 2>&1
