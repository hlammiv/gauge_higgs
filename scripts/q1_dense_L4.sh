# Dense q=1 search on L=4 (fast/exploratory) USING THE MULTI-TIMESCALE INTEGRATOR.
# Finer grid than the L=8 pilot-plus (kappa step 0.05, wide beta) to resolve the
# boundary ridges; n_scalar=6 (Sexton-Weingarten) + auto-tune keep every point cheap
# and stable (the stiff deep-Higgs points run at low nmd instead of nmd~40).
# CAVEAT: L=4 in D=4 has large finite-size effects -- the Coulomb phase / beta_c are
# quantitatively distorted (beta_c shifts from ~1.01); this is a QUICK exploratory map
# of the structure, not a quantitative result (use L=8/12/16 + FSS for that).
LS=(4)
QS=(1)
BETAS=($(seq -f "%.2f" 0.80 0.05 1.30))   # 11 values
KAPPAS=($(seq -f "%.2f" 0.00 0.05 0.60))  # 13 values
LAMBDA=0.5
NTHERM=400
NMEAS=1500
NMD=6              # starting nmd; auto-tune + n_scalar keep it low
TAU=1.0
MEASURE_EVERY=2
AUTOTUNE=1
N_SCALAR=6         # multi-timescale: 6 scalar sub-steps per gauge step
CAMPAIGN_SEED=20260601
NDIM=4
