# Clean q=1 "pilot-plus" config for scripts/u1_scan_campaign.sh
# A trustworthy single q=1 (beta,kappa) phase diagram: auto-tuned nmd (stable
# everywhere), monopole column on (sharp chi_rho Coulomb-confinement locator),
# link susceptibility (Higgs-confinement), finer 7x7 grid focused on beta_c~1.0-1.05,
# and ~10x the pilot's statistics. Validates that the boundaries sharpen as expected
# before committing to the full q-series (#28).
LS=(8)
QS=(1)
BETAS=(0.90 0.95 1.00 1.05 1.10 1.15 1.20)
KAPPAS=(0.0 0.10 0.20 0.30 0.40 0.50 0.60)
LAMBDA=0.5
NTHERM=1000
NMEAS=4000
NMD=10            # starting nmd; auto-tune adjusts per point
TAU=1.0
MEASURE_EVERY=2
CAMPAIGN_SEED=20260530
NDIM=4
