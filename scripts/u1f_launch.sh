#!/usr/bin/env bash
# Persistent core-filling launcher for a frozen-U(1) campaign queue.
# usage: scripts/u1f_launch.sh <queue_file> <ncores> [outdir]
# Detaches via setsid so it survives an ssh disconnect; resumable (u1f_job.sh skips finished jobs).
set -u
Q="$1"; N="$2"; OUT="${3:-u1f_campaign}"
mkdir -p "$OUT"
U1F_OUT="$OUT" setsid bash -c "xargs -P $N -L1 -a '$Q' bash scripts/u1f_job.sh" \
  > "$OUT/launch_$(basename "$Q").log" 2>&1 < /dev/null &
sleep 1
echo "launched: $(wc -l < "$Q") jobs, $N cores, out=$OUT"
