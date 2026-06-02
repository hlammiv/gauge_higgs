# Stage B -- m_gamma PILOT on the boundary charges, to RIGOROUSLY test the q>=5
# Coulomb wedge (the headline). Full (beta,kappa) plane to kappa=1.0 WITH m_gamma
# measured, so the Coulomb phase (massless photon) is identified definitively rather
# than inferred. Charges chosen to bracket the physics:
#   q=1  -- light / Fradkin-Shenker (Higgs-confined connected; sanity check)
#   q=4  -- last charge expected to have a FINITE triple point (Z_4 deep-Higgs)
#   q=5  -- first charge expected to RECEDE (Z_5 deep-Higgs has a Coulomb phase)
#   q=8  -- deep in the receding regime
# If q=4 shows a finite triple point and q=5,8 show m_gamma staying small to kappa=1.0
# (Coulomb wedge persists), that's the q<=4-finite / q>=5-receding result.
LS=(8)
QS=(1 4 5 8)
BETAS=($(seq -f "%.2f" 0.80 0.10 1.30))    # 6
KAPPAS=($(seq -f "%.2f" 0.00 0.10 1.00))   # 11
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
CAMPAIGN_SEED=20260530   # distinct from scenarioA -> independent streams
NDIM=4
# 4 q x 6 beta x 11 kappa = 264 points WITH m_gamma (cor_*.dat + m_gamma columns).
