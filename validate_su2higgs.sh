#!/usr/bin/env bash
# Finite-size validation of the SU(2)-fundamental Higgs confinement<->Higgs transition.
# Scans kappa across the transition on two volumes and records the gauge-invariant link
# energy L_link, its susceptibility chi_link = V*Var(L_link), Binder cumulant, and the
# plaquette. The susceptibility peak locates the transition and should grow/sharpen with
# volume (finite-size scaling). beta=2.3, lambda=0.5 fixed.
set -u
cd "$(dirname "$0")"
BIN=./build/hmc_higgs
OUT=validation_su2higgs.dat
THREADS=${OMP_NUM_THREADS:-16}
NTHERM=120; NMEAS=300; NMD=18; SEED=1234
printf "# SU(2) fund Higgs  beta=2.3 lambda=0.5  ntherm=%d nmeas=%d nmd=%d\n" $NTHERM $NMEAS $NMD > $OUT
printf "# %-3s %-7s %-9s %-9s %-9s %-9s %-7s %-6s\n" L kappa L_link err chi_link chi_phi Binder acc >> $OUT
for Lx in 6 8; do
  for k in 0.18 0.21 0.23 0.245 0.255 0.265 0.28 0.31; do
    o=$(OMP_NUM_THREADS=$THREADS $BIN fund $Lx 2.3 $k 0.5 $NTHERM $NMEAS $NMD 1.0 $SEED)
    ll=$(  echo "$o" | awk '/^L_link/   {print $3}')
    er=$(  echo "$o" | awk '/^L_link/   {print $5}')
    cl=$(  echo "$o" | awk '/^chi_link/ {print $3}')
    cp=$(  echo "$o" | awk '/^chi_phi/  {print $3}')
    bn=$(  echo "$o" | awk '/^Binder/   {print $3}')
    ac=$(  echo "$o" | awk '/^acceptance/{print $3}')
    printf "%-5s %-7s %-9s %-9s %-9s %-9s %-7s %-6s\n" "$Lx" "$k" "$ll" "$er" "$cl" "$cp" "$bn" "$ac" | tee -a $OUT
  done
done
echo "=== DONE ==="
