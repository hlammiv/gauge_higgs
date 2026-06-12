#!/usr/bin/env bash
# V(R) WEDGE-COLUMN SWEEP (frozen `pot` mode) -- the position-space static-potential cross-check of the
# structure-factor m_gamma wedge. At each (q,beta,kappa) it emits the full charge-1 Wilson grid W[R][T];
# scripts/u1f_vr.py then reconstructs V(R) and Cornell-fits c+sigma*R-alpha/R. The SHAPE labels the phase:
#   Confined sigma>0 (linear) | Coulomb sigma~0,alpha>0 (-1/R tail, UNSCREENED, wedge) | Higgs flat (SCREENED).
#
# Defaults = the near-wall wedge column where V(R) is clean (Higgs screening length < 1 lattice spacing):
#   q in {2,4,5,6,8} x beta in {1.2,1.5} x kappa in {1.0,1.5,2.5}  = 30 jobs, L=16, Rmax=8, nmeas~625.
# Single-threaded per job (sweep is serial), parallelised across cores -- ~10-15 min/job, ~30 jobs / 18 cores
# ~ 2 waves ~ 30-40 min on lucia. Detaches via setsid (survives disconnect); resumable (skips finished jobs).
#
# usage:  scripts/u1f_vr_sweep.sh                          # launch on lucia (18 cores -> u1f_vr/)
#         OUT=u1f_vr NCORES=18 scripts/u1f_vr_sweep.sh
#         KAPPAS="0.6 1.0 1.5 2.0 2.5" scripts/u1f_vr_sweep.sh   # finer kappa column
# analyze: python3 scripts/u1f_vr.py u1f_vr/*.out
set -u
OUT="${OUT:-u1f_vr}"; NCORES="${NCORES:-18}"
L="${L:-16}"; RMAX="${RMAX:-8}"
NSWEEP="${NSWEEP:-2500}"; NTHERM="${NTHERM:-1000}"; MEAS="${MEAS:-4}"; NOR="${NOR:-3}"
QS="${QS:-2 4 5 6 8}"; BETAS="${BETAS:-1.2 1.5}"; KAPPAS="${KAPPAS:-1.0 1.5 2.5}"
BIN="${U1F_BIN:-./build/u1_frozen}"
mkdir -p "$OUT"

if [ ! -x "$BIN" ]; then echo "ERROR: binary $BIN not found/executable (build it first)"; exit 1; fi

# build the queue (one `pot` CLI per line); deterministic seed per point so a resume maps to the same file
Q="$OUT/queue.txt"; : > "$Q"
i=0
for q in $QS; do for b in $BETAS; do for k in $KAPPAS; do
  seed=$((70001 + i)); i=$((i + 1))
  echo "pot $L $b $k $q $RMAX $NSWEEP $NTHERM $MEAS $seed $NOR" >> "$Q"
done; done; done
N=$(wc -l < "$Q")
echo "queue: $N pot jobs  L=$L Rmax=$RMAX nsweep=$NSWEEP nmeas~$((NSWEEP / MEAS))  ->  $OUT/  on $NCORES cores"
sed 's/^/   /' "$Q"

# launch detached; u1f_pot_job.sh skips already-finished outputs
POT_OUT="$OUT" U1F_BIN="$BIN" setsid bash -c \
  "xargs -P $NCORES -L1 -a '$Q' bash scripts/u1f_pot_job.sh" \
  > "$OUT/launch.log" 2>&1 < /dev/null &
sleep 1
echo "launched (detached, $NCORES cores)."
echo "  progress:  echo \"\$(grep -lc '# u1_frozen pot:' $OUT/job_*.out 2>/dev/null | wc -l) / $N done\""
echo "  analyze :  python3 scripts/u1f_vr.py $OUT/job_*.out"
