# Scenario A -- L=8 cartography for the full charge series q=1..8 (D=4).
#
# SPARSE (beta,kappa) grid by design (NOT dense): it brackets the Coulomb<->confined
# gauge line (beta_c ~ 1.01, ~q-independent) and the Higgs<->confined line / triple
# point; REWEIGHTING (reweight.hpp) then interpolates between nodes. Dense grids are
# the wrong tool for beta_t/kappa_t precision -- the right tool is sparse grid + high
# stats + reweighting (+ L-ladder FSS as a separate follow-up).
#
# Deliverable: beta_t(q) for q=1..8 and the kappa_t(q) trend (finite & growing for
# q<=4; Coulomb wedge receding to kappa->inf for q>=5) at ONE volume (L=8). The
# quantitative triple-point localization that needs continuum/FV control is the
# targeted L={12,16} FSS stage, run afterward.
#
# Multi-timescale (n_scalar=6) + per-(L,q) PRETUNE keep every point cheap and stable
# (nmd stays ~low and flat across the plane; no per-point autotune calibration blowup).
#
# Run (lenore, 32 cores -> JOBS=16 auto):
#   scripts/u1_scan_campaign.sh --config scripts/scenarioA.sh --out u1scan_scenarioA          # DRY-RUN
#   scripts/u1_scan_campaign.sh --config scripts/scenarioA.sh --out u1scan_scenarioA --run    # EXECUTE
LS=(8)
QS=(1 2 3 4 5 6 7 8)
BETAS=($(seq -f "%.2f" 0.80 0.10 1.30))    # 6 values: 0.80 .. 1.30 (brackets beta_c~1.01)
KAPPAS=($(seq -f "%.2f" 0.00 0.10 0.60))   # 7 values: 0.00 .. 0.60 (kappa=0 = pure-gauge reference)
LAMBDA=0.5
NTHERM=500
NMEAS=2000
NMD=6                # pretune START value (final nmd set per (L,q) by pretune)
TAU=1.0
MEASURE_EVERY=2      # -> 4500 traj/point (500 therm + 2000 x 2)
AUTOTUNE=0           # ignored under PRETUNE=1; explicit for clarity
N_SCALAR=6           # multi-timescale Sexton-Weingarten scalar sub-steps/gauge step
PRETUNE=1            # tune nmd once per (L,q) at the stiff corner, reuse across the grid
PRETUNE_NTHERM=150
CAMPAIGN_SEED=20260529
NDIM=4
# 8 q x (6 beta x 7 kappa) = 336 points x 4500 traj ; ~4 h wall on lenore (16 jobs).
