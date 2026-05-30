# Stage-0 VALIDATION PILOT config for scripts/u1_scan_campaign.sh
# Goal: end-to-end pipeline validation + preliminary signals, NOT production.
#   - kappa=0 column  -> pure compact U(1): locate the bulk beta_c ~ 1.0-1.01
#     (plaquette / monopole-density jump; first with-HMC test of rho_M, m_gamma, Creutz).
#   - q=1, (beta,kappa) ~ (0.80-0.90, 0.45-0.60): bracket the Fradkin-Shenker
#     endpoint D ~ (0.8485, 0.526) (link-energy susceptibility).
# Coarse single L=8 grid; scale up L (add 12,16), NTHERM/NMEAS, and refine the grid
# for production once the pipeline + signals look right.
LS=(8)
QS=(1)
BETAS=(0.80 0.90 0.95 1.00 1.05 1.10)
KAPPAS=(0.0 0.30 0.45 0.53 0.60)
LAMBDA=0.5
NTHERM=400
NMEAS=800
NMD=20
TAU=1.0
MEASURE_EVERY=2
CAMPAIGN_SEED=20260529
NDIM=4
