# What We Give Up By Freezing the Radial Mode (|phi|=1, London limit lambda->inf)

**Question:** We want to freeze the Higgs modulus (|phi|=1, the pure-phase / London / lambda->inf limit) to
escape the deep-Higgs first-order matter sampling wall. What physics do we actually give up, and is the prior
claim that the triple point's "EXISTENCE, ORDER, and q-dependence are invariant under freezing" true?

**Short answer:** Freezing removes exactly one coherent thing — the dynamical Higgs **amplitude** mode — and
every loss is a corollary of that. **None** of the losses touch the #28 triple-point locators, the #26 q>=5
Coulomb-wedge locators, or the kappa->inf Z_q digitization headline; those are all gauge/phase-sector facts that
the frozen model reproduces faithfully (and where the FS/Bowler theorems were literally proven *in* the frozen
model). The "EXISTENCE and q-dependence invariant" half of the prior claim is **correct**. The "**ORDER**
invariant" half is **overclaimed**: the order of the Coulomb-Higgs transition is lambda-tuned through a
**tricritical point** and freezing pins us on one side of it. It is *safe at our operating point* (lambda=0.5),
but false as a blanket statement. The genuine casualties are off-path tasks: **#16 Higgs spectroscopy** and
**#25 type-I/II superconductor** physics, which are frozen-incompatible by construction and must run on the
dynamical code.

---

## 1. LEDGER — everything given up, sorted so what matters is at the top

Legend: "hurts #28/#26?" = does it blind a triple-point / Coulomb-wedge **locator** or change the existence/
q-dependence claim. "hurts digitization?" = does it affect the kappa->inf Z_q endpoint (group, beta_c, order).

### TIER A — MATTERS / must be qualified (read these)

| # | What is given up | Hurts #28/#26? | Hurts digitization? | Recoverable how? |
|---|---|---|---|---|
| **A1** | **The blanket "ORDER invariant" claim itself.** The order of the Coulomb-Higgs transition is lambda-dependent through a **tricritical point** at lambda_t ~ 0.005 (Espriu/Aguado et al., hep-lat/9708011, 24^4): first order below, second order above. Frozen = lambda=inf = locked on the **second-order** side. A type-I (small-lambda) physical target is first order while its frozen surrogate is second order -> **order flips**. lambda_t grows toward *lower* beta (nearer the triple point), so the disagreement window is *widest* exactly where #28 lives. | **PARTIAL** — does NOT move the triple-point existence/intersection, but if the paper quotes the *order* of the Higgs/Coulomb arm, the frozen number is the lambda=inf statement, not the lambda=0.5 one. SAFE at our lambda=0.5 (~100x above lambda_t, both frozen and dynamical on the 2nd-order side). | **NO** — the kappa=inf Z_q order (1st-order q<=4, two transitions q>=5) is lambda-independent by the unitary-gauge argument. | **One dynamical lambda=0.5 control slice** cutting the Coulomb-Higgs line near the triple point: Binder cumulant / Delta-F FSS, frozen vs dynamical, on a single line. This is the single cheapest insurance run for the paper. **Caveat:** do NOT read order from L<=8 hysteresis — hep-lat/9708011 shows hysteresis *misreads* order in the thermodynamic limit for BOTH models; need Binder/Delta-F at L>=16. |
| **A2** | **The Higgs (amplitude) boson and its mass m_H** (task **#16** spectroscopy). The radial fluctuation rho-<rho> *is* the physical Higgs; m_H^2 ~ 8 lambda <rho>^2 at tree level. Freezing sends it to the cutoff: the frozen scalar is a pure compact phase angle, so every |phi|-amplitude correlator is identically trivial. The 0^++ scalar channel of the H-multiplet is gone. | **NO** — no #28/#26 locator is the Higgs amplitude pole. | **NO** — at kappa->inf the radial mode decouples anyway (effective mass grows ~D*kappa); the Z_q endpoint never contained it. | **Dynamical-modulus only.** Make #16 an explicitly **dynamical** task: radial timeslice correlator C(dt) of (|phi_x|-<|phi|>) at a few (beta,kappa) in the Higgs phase at *moderate* kappa where HMC is healthy. Cannot be obtained from any frozen run. |
| **A3** | **The Ginzburg-Landau ratio m_H/m_W and the type-I vs type-II classification** (task **#25**). m_W = sqrt(2 kappa) q SURVIVES (gauge sector, eaten-Goldstone mass — gh_vecmass/gaugeboson_op). m_H does not. So frozen pins kappa_GL = m_H/(sqrt2 m_W) -> **infinity** = extreme type-II / London corner ONLY. The type-I regime and the type-I<->type-II crossover line are **unreachable**. #25 in its full superconductor form (B_c1/B_c2, resolved vortex cores, Meissner-vs-vortex) is **frozen-incompatible by construction**. | **NO** — type-I/II is orthogonal to the triple point (topology/screening, not kappa_GL). | **NO** — kappa_GL classifies the broken Higgs phase; the discrete Z_q residual has no amplitude mode. | **PARTIAL / dynamical-only for the NUMBER.** m_H/m_W cannot come from any frozen run (m_H=inf). Run #25 as its own targeted **finite-lambda** dynamical campaign: vary lambda (small ~0.1 -> type-I/first-order radial; ~O(1) -> type-II) at moderate kappa. Frozen runs *can* still do the type-II side (vortex lattices, flux quantization in 2pi/q units, Meissner expulsion = penetration depth = gauge sector). Flag #25 "frozen-incompatible by construction." |

### TIER B — real losses, but provably off the #28/#26/digitization path

| # | What is given up | Hurts #28/#26? | Hurts digitization? | Recoverable how? |
|---|---|---|---|---|
| **B1** | **The small-lambda (lambda<=0.13) fluctuation-induced FIRST-ORDER radial line** (Damgaard-Heller; the HLM/Coleman-Weinberg mechanism). At beta=0 the gauge field integrates out exactly and for lambda<=0.13 there is a first-order jump in <rho^2>, **Q-independent**, driven by non-trivial minima of V_eff(rho)=lambda(rho^2-1)^2+rho^2-ln(rho)-(1/4)ln I_0(2 kappa rho^2); it "persists for all beta" and can join the Coulomb-Higgs line, *cutting the confined phase in two*. This entire boundary is **identically absent** frozen — D-H's own "first instance of a difference to the fixed-length limit." | **PARTIAL** — does NOT exist at our default lambda=0.5 (0.5 >> 0.13), so #28/#26 at lambda=0.5 lose nothing. It only forbids claiming the frozen map represents the small-lambda corner of the (beta,kappa,lambda) cube. NB: this radial line *is* the deep-Higgs B-jump that defeated HMC/LLR/muca — freezing escapes the wall by *deleting* it, which is legitimate for the diagram but means we never map it. | **NO** — it lives at small kappa and never reaches the kappa=inf boundary (D-H Fig.3: deep-kappa Z_2 line hits 0.4407 even at lambda=0.01). | **Dedicated small-lambda dynamical scan** (lambda=0.01-0.13, beta~0), with D-H's analytic V_eff(rho) as the oracle (matches MC "almost perfectly"). Cheap; not on the current task list. |
| **B2** | **The coherence length xi = 1/m_H** (-> 0 in the London limit) and **resolved vortex-core structure**. A frozen vortex has phase winding and flux but NO amplitude core; D-H's vortex signature (rho-dip exactly where local flux peaks) is gone, so xi and the xi-vs-penetration-depth length-scale separation are unmeasurable. | **NO** — flux/screening locators (sigma_1, m_gamma) sense the phase winding + penetration depth (gauge sector, survives), not the amplitude core. | **NO** — xi->0 is the controlled idealization OF the limit we digitize; a Z_q theory has no coherence length. | **Rides on the same dynamical slice as m_H** (xi=1/m_H). One radial-correlator measurement gives both. |
| **B3** | **The external-field superconductor program** (Damgaard-Heller Sec.3, Figs.5-10): B_c1/B_c2, Meissner vs vortex phase, vortex identification by cooling. Two ingredients: (i) the external-field plaquette shift theta^ext = 2pi n_ex/L^2 is a gauge add-on, untouched by freezing (not currently coded — grep shows only monopole/photon flux, no n_ex); (ii) the type-I/II identification itself needs the dynamical modulus (rho-dip core, finite kappa_GL). | **NO** — separate external-field program, orthogonal to the zero-field (beta,kappa) diagram #28/#26 live on. | **NO** — not part of the kappa->inf reduction. | **PARTIAL.** Meissner expulsion + flux quantization: yes even frozen (gauge sector + ~tens of lines to add theta^ext to u1.hpp). Type-I-vs-II, B_c1/B_c2 separation, resolved cores: **NO** — need the dynamical modulus. This is the largest single program freezing forecloses, but it is a deliberately-separate sub-project = #25 in full form. |

### TIER C — bookkeeping / map distortion (no claim at risk, but label everything)

| # | What is given up | Hurts #28/#26? | Hurts digitization? | Recoverable how? |
|---|---|---|---|---|
| **C1** | **Quantitative coordinates.** kappa_frozen ~ kappa*<rho^2> is a *non-uniform, discontinuous* reparametrization, not a global rescale: <rho^2>(beta,kappa,lambda) varies across the diagram and JUMPS at first-order lines. At lambda=0.5 the classical minimum sits at rho^2 = 1-1/(2lambda) = 0, so <rho^2> is entirely hopping-generated (~1.9 at kappa=0.3, ~3.8 at 0.5, ~8 at kappa=1). Frozen (beta_t,kappa_t) numbers are NOT the lambda=0.5 numbers; line shapes/intersection angles deform. | **PARTIAL** — existence/q-dependence/topology are map-invariant; every quoted coordinate must be labeled "frozen (London-limit) model" (standard 1980s practice; Bowler and D-H refs [5,6] all publish frozen coordinates). | **NO** for the claim; YES for bookkeeping — the kappa->inf asymptote (group, beta_c, order) is identical, only the finite-kappa trajectory is reparametrized. | Measure <rho^2>(beta,kappa) on one dynamical lambda=0.5 slice; overlay frozen diagram under kappa -> kappa<rho^2>. D-H eq.(2) beta_Q = beta(1 - c/(kappa rho_bar^2 Q^2)) predicts the leading deformation analytically. |
| **C2** | **Rate of approach to the Z_q limit.** Dynamical (rho_bar^2 ~ D kappa/lambda grows) converges ~1/kappa^2; frozen (rho_bar^2=1) converges ~1/kappa (Bowler eq.5: beta_eff = beta_p(1 - 1/(2 beta_L)) for q=2). So the frozen model needs LARGER raw kappa for the "large-kappa slice == pure Z_q" matching test. | **NO** — a compute-budget statement; the frozen sampler (heatbath+overrelax+Z_q link jump) is exact at any kappa, so large kappa is cheap there, unlike HMC. | **PARTIAL, cuts both ways.** Cost: slower raw-kappa convergence. **Benefit:** the approach becomes a clean one-parameter 1/kappa series with NO lambda-dependent prefactor (Bowler eq.5, literature-validated; q=2 demonstrably -> Z_2 self-dual 0.4407), making the sigma_1/|P_1| matching extrapolation *cleaner*, not weaker. | Quote matching as sigma_1(beta,kappa) -> sigma_1^{Zq}(beta_q(kappa)) with beta_q from the frozen 1/kappa series; cross-check one point on the dynamical code at moderate kappa. |
| **C3** | **rho=0 amplitude zeros and latent-heat NUMBERS.** rho is a **gauge singlet** (U(1) acts only on the phase), so freezing it **cannot** change WHICH discrete group is residual (= kernel of charge-q rep = Z_q, pure representation theory). The only loophole is rho_x=0 sites (stabilizer jumps to U(1)); their weight vanishes at kappa->inf (each costs ~kappa hopping). Frozen latent heats / B-jump magnitudes differ from lambda=0.5 — but the dynamical B-jump (high branch ~6000 vs frozen ceiling 2*vol*D) is dominated by radial self-amplification *artifact*, not Z_q physics. | **NO** for existence; the deep-kappa first-order latent-heat numbers quoted frozen (bounded B) are the well-posed object. | **NO — strengthens it.** "Does freezing change the residual group?" is rigorously **NO**; freezing even *removes* the only configs (rho=0 cores) that formally muddy the unitary-gauge argument at finite kappa. | Vortex-core/amplitude-zero physics belongs to the dynamical #25 program at moderate kappa, if ever wanted. |

---

## 2. Verdict on "ORDER is invariant under freezing"

**REFUTED as a blanket claim; SAFE at our operating point.** It must be downgraded, not accepted.

**Where it HOLDS (rigorously):**
- **At the digitization endpoint kappa->inf:** pure Z_q gauge theory, lambda-independent. Damgaard-Heller state the
  unitary-gauge argument *in the variable-length model* (links lock to e^{2pi i n/q} for ANY lambda); Fradkin-Shenker
  prove the same in the fixed-length model. Order there is a Z_q fact: single first-order for q<=4, two transitions
  (intermediate massless phase) for q>=5. The kappa=inf beta_c is lambda-invariant across lambda=0.01,1.0,inf
  (D-H Fig.3: deep-kappa Z_2 line hits (1/2)ln(1+sqrt2)=0.4407 even at lambda=0.01; Bowler q=2 -> same 0.4407).
  **This carries #26 and the digitization headline cleanly.**
- **On the deep-kappa Z_q continuation lines:** order inherited from the endpoint via beta_q(kappa); Bowler Fig.4
  shows the STRONG first-order at (beta_p=0.5, beta_L=4.0) *in the frozen model*.
- **At our actual operating point lambda=0.5:** ~100x above the tricritical lambda_t ~ 0.005. lambda=0.5 (dynamical)
  and lambda=inf (frozen) are BOTH on the second-order side, so for the runs actually planned they agree.
- **EXISTENCE / connectivity / q-dependence:** FS is itself a fixed-length theorem; the residual group is pure
  representation theory (gauge-singlet rho). These are frozen-robust. Note also: the frozen 4D model is NOT
  order-trivial — gauge-fluctuation-induced (Coleman-Weinberg) weak first order *survives* freezing (FS Fig.2:
  the q=1 XF line is first-order in the fixed-length model). The naive "frozen = always second order" is itself wrong.

**Where it FAILS:**
- **The Coulomb-Higgs transition has a tricritical point in lambda** at lambda_t ~ 0.005 (hep-lat/9708011, 24^4):
  first order below, second order above. Freezing (lambda=inf) sits on the second-order side; a **type-I target
  (lambda < lambda_t) is first order while its frozen surrogate is second order -> the order FLIPS.** And lambda_t
  *increases* toward the lower beta of the triple point, widening the disagreement window exactly at #28.
- **Damgaard-Heller's lambda<=0.13 first-order radial line** exists dynamically (joins the Coulomb-Higgs line at
  small beta) and is **identically absent** frozen -> even **EXISTENCE** of a phase boundary is not invariant in
  the small-lambda corner. D-H: "differs significantly from the fixed-length limit."
- **Mechanism** = exactly the type-I/II / Halperin-Lubensky-Ma lever: integrating out gauge fluctuations generates
  a |phi|^3 term driving a fluctuation-induced first-order transition in type-I (small-lambda, kappa_GL<1/sqrt2),
  which weakens to continuous/inverted-XY in deep type-II (London). Freezing pins kappa_GL -> inf (extreme type-II),
  discarding the entire first-order type-I half of the line.

**Corrected statement to put in the paper:**
> "q-dependence and the kappa=inf Z_q endpoint order are lambda-invariant. EXISTENCE/topology and the order of the
> Coulomb-Higgs arm are invariant only for lambda above the (beta-dependent) tricritical threshold ~0.005 — i.e. at
> our lambda=0.5, NOT in the type-I / small-lambda corner."

**Methodological caveat that bites in BOTH models:** hep-lat/9708011 prove small-lattice double-peak/hysteresis
signals vanish in the thermodynamic limit down to lambda~0.005, so small-lattice "first order" reads were wrong.
The frozen-radial-route plan to read order by L=6-8 hysteresis/thermal-cycling is exactly the method that paper
refutes. **Determine the Coulomb-Higgs order by Binder/Delta-F FSS at L>=16, never by L<=8 hysteresis** — this is a
problem with the *method*, independent of freezing.

---

## 3. Bottom line

**FREEZE — for the #28/#26 cartography and the kappa->inf Z_q digitization headline.** The decision is rigorously
justified: the kappa=inf endpoint group, its beta_c, its order (q<=4 single first-order / q>=5 two-transition split),
and the full q-dependence are all lambda-independent (D-H unitary-gauge argument in the variable-length model; FS in
the fixed-length model; Bowler/Fernandez-Sudupe validations), and the residual group cannot change because rho is a
gauge singlet. Every #28/#26 *locator* (sigma_1, |P_1|, m_gamma, chi_link) is gauge/phase-sector and frozen-faithful.
What freezing genuinely deletes — the Higgs amplitude m_H (#16), the GL ratio m_H/m_W and type-I/II classification
(#25, frozen-incompatible by construction since m_H=inf), xi=1/m_H, vortex cores, and the small-lambda HLM first-order
radial line — is all **off the #28/#26/digitization path** and recoverable only on the dynamical code, which we keep.

**Two things the lead must do:** (1) **Drop or qualify the word "ORDER"** in the invariance claim — it holds for the
endpoint/q-dependence and at lambda>=0.5, but NOT in the type-I corner. (2) Run **#16 and #25 as a separate dynamical-
modulus, finite-(and small-)lambda sub-project** — they are physically lost by freezing, not recoverable from any
frozen run.

**The ONE control measurement that buys back the most:** **one dynamical-modulus lambda=0.5 control slice** that
(a) measures <rho^2>(beta,kappa) to pin the kappa_frozen ~ kappa*<rho^2> map, (b) extracts m_H and m_H/m_W (gives xi
for free), and (c) cross-checks the Coulomb-Higgs **order** frozen-vs-dynamical near the triple point via **Binder/
Delta-F FSS at L>=16** (not hysteresis). That single slice converts the contested "order invariant" assertion into a
measured cross-check and closes the referee's strongest objection without surrendering the frozen-route speedup for
the production diagram.

**Readiness flag:** the U(1) frozen path is **not yet built**. Frozen-length is implemented only in the SU(N)/general
driver (src/hmc/gauge_higgs_hmc.hpp geodesic_drift + project_pi_tangent, GH_FROZEN env in src/hmc_higgs_multi.cpp);
**src/u1_scan.cpp has no frozen mode**. Any prior statement that "frozen-length is implemented" overstates U(1)
readiness — a frozen U(1) heatbath/overrelax (+ Z_q link jump) sampler still needs to be written before the frozen
production campaign can run.
