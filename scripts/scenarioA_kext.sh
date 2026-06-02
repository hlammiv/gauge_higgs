# Scenario A -- EXTENDED kappa to 1.0. Resumes into u1scan_scenarioA (skips the
# kappa<=0.6 points already done) and ADDS kappa in {0.7,0.8,0.9,1.0} so that:
#   - every charge's Higgs<->confined line is bracketed (q=2's sat at the old 0.6 edge);
#   - the q>=5 Coulomb wedge becomes visible -- does the Coulomb phase persist up to
#     large kappa (no triple point, the headline) or pinch off?
# New points: 8 q x 6 beta x 4 new kappa = 192 (the existing 336 are reused as-is).
# Pre-tune now tunes at the stiffer corner (kappa=1.0); high-q/high-kappa points will
# get larger nmd -- check u1scan_scenarioA/pretune.log before trusting the wall estimate.
LS=(8)
QS=(1 2 3 4 5 6 7 8)
BETAS=($(seq -f "%.2f" 0.80 0.10 1.30))    # 6 values: 0.80 .. 1.30
KAPPAS=($(seq -f "%.2f" 0.00 0.10 1.00))   # 11 values: 0.00 .. 1.00 (extended from 0.60)
LAMBDA=0.5
NTHERM=500
NMEAS=2000
NMD=6
TAU=1.0
MEASURE_EVERY=2
AUTOTUNE=0
N_SCALAR=6
PRETUNE=1
PRETUNE_NTHERM=150
CAMPAIGN_SEED=20260529   # same as scenarioA -> deterministic; resume keys on (beta,kappa) dirs
NDIM=4
