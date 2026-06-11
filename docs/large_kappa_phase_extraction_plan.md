# Large-κ Gauge+Higgs Phase Extraction — Decision-Ready Roadmap

**Scope:** U(1)+charge-q Higgs (β = gauge coupling, κ = Higgs hopping). Settle (#28) the
finite-κ triple point for q ≤ 4, (#26) the q ≥ 5 intermediate-Coulomb wedge, and the
κ→∞ reduction to the discrete residual Z_q gauge theory.
**Status:** PLANNING ONLY — no production runs without the lead's explicit go-ahead.
**Author:** technical-lead integration of 5 candidate designs + their 3-lens adversarial verdicts.

---

## 1. Executive summary

The single approach most likely to actually crack the large-κ matter axis is the **dual /
worldline (flux) representation** of the U(1)+charge-q Higgs model: an exact character-expansion
rewrite of *our* cosine action (Bessel `I_n(β)` plaquette weights, integer link fluxes `k`, dimers
`l`, positive radial site weights `W(f)`) in which the deep-Higgs reaction coordinate becomes an
**integer flux total `F` with ⟨B⟩ = ⟨F⟩/κ exact**, sampled by O(1)-cost integer Metropolis with a
multicanonical weight in the **accept ratio** — so all three documented failure modes (HMC can't
tunnel, LLR can't hold high B, muca-in-force froze at κ_eff≈−134) are absent *by construction*, not
by tuning. It was the only design all three verifier lenses scored ≥7 and the only one a verifier
called "the first of the five methods that should actually work." Its honest gaps — the gauge
A-sector still needs its own integer muca, the σ₁/|P₁| defect estimator at high q is the
least-certain payoff, and L=16 order-FSS round-trip cost must be measured — are additive fixes, not
architectural blockers.

The **fallback / parallel-insurance track**, recommended to start *concurrently* because it is far
cheaper and de-risks the dual, is **accept-side multicanonical-in-B** (g(B) in the Metropolis accept,
never the force) — the explicit recorded fix that has never been tried; it shares the same MucaB
container and the same `scripts/u1_llr_maxwell.py` analysis, and is the natural sampler for the
*moderate*-κ triple-point neighborhood where the barrier is weakest. The **third, near-zero-risk
track that we build regardless** is the κ→∞ matching infrastructure: a standalone **pure-Z_q gauge
driver** plus an **exact dB=0 Z_q-hop link move** (`θ→θ+2πm/q`, the abelian GH_GUPD generalization)
that restores gauge ergodicity on large-κ matching slices without ever touching the B-jump. The key
decisions: **(1)** commit to the dual as primary but gate its L=16/high-q claims on measured
round-trip tests; **(2)** run accept-side muca as cheap insurance + the moderate-κ triple-point tool;
**(3)** build the pure-Z_q driver + hop move now — they are required for the headline κ→∞ figure and
carry no crux risk.

---

## 2. The problem, precisely

**Action (authoritative, identical across `u1.hpp`/`scan_obs.hpp`/`u1_llr.hpp`):**
`S = β·A − κ·B + S_pot`, sampled `exp(−S)`, with
- `A = Σ_plaq (1 − cos θ_plaq)` — extensive, conjugate to β.
- `B = Σ_{x,μ} 2 Re[conj(φ_x) e^{iqθ} φ_{x+μ}]` — extensive hopping sum; the conjugate variable is **−B**.
- `S_pot = Σ_x [|φ|² + λ(|φ|²−1)²]`.

Reweighting coords `(E1,E2)=(A,−B)`, `(λ1,λ2)=(β,κ)`.

**The gauge axis is solved.** The pure-gauge compact-U(1) bulk transition sits at β_c ≈ 0.9–1.0
(literature ~1.01); LLR-in-A crosses its first-order barrier with ⟨A⟩(β) agreeing with direct HMC to
mean|dA| = 0.5–1.1%. χ_plaq / σ_fund / Creutz ratios / |P_fund| are all wired and validated.

**The matter axis is the open crux.** The deep-Higgs (large-κ) symmetric↔Higgs transition is a
strongly first-order jump in B whose latent-heat × volume barrier has defeated **four** distinct
samplers, established rigorously and not to be relitigated:
1. **Local HMC** cannot tunnel it (cost ∝ latent-heat × V); acceptance stays 0.96–1.0, so it is a
   sampling limit, not an integrator error. Hot−cold ⟨plaq⟩ gap grows with q (q=1→0.00 control,
   q=4→0.48, q=8→0.62).
2. **2D-LLR matter windowing** cannot HOLD high B: guarded-HMC acceptance collapses (<0.3),
   res_E2 ≈ 3× window, 43/50 high-B cells hit the nwiden=3 ceiling. Budget-independent.
3. **Naive muca-in-FORCE** (`keff()=κ+g'(B)` in the HMC kick, `u1.hpp:201`) froze: WL g ran away
   (~6700/bin), κ_eff ≈ −134, HMC rejected everything.
4. **Single-axis β/κ replica tempering**: swaps die mid-ladder at the extensive-variance wall
   (killer pairs (2,4)=(4,6)=0.0000 swap acceptance).

Both headline questions (#28 finite-κ triple point; #26 q≥5 Coulomb wedge) and the κ→∞ discrete
reduction live in this broken matter direction. The certificate also needs latent jumps ΔA, ΔB for
the Clausius–Clapeyron check `dκ/dβ = −ΔA/ΔB`, where ΔB is precisely the unsampleable quantity.

A second, separable hole: **ρ_M (monopole density) and m_γ (photon mass) are BLIND** to the
deep-Higgs Z_q-confined phase (charge-q condensate Higgs-screens monopoles; m_γ is massive in both
Higgs and Z_q-confined). The unscreened locators are the **charge-1 string tension σ₁** and
**charge-1 Polyakov |P₁|** (q=1 = the screened negative control).

---

## 3. Method assessment table

Scores are the mean of the three adversarial-lens verdicts (Tunneling-correctness / Implementation /
Physics-payoff). "Will-it-tunnel" = does it cross the deep-Higgs B-jump.

| # | Method | Will it tunnel the B-jump? | Build effort | Physics payoff | Lens scores | Mean | Verdict |
|---|--------|---------------------------|--------------|----------------|-------------|------|---------|
| 1 | **Dual / worldline flux rep** | **YES** — F is integer, ⟨B⟩=⟨F⟩/κ exact, free q-independent l-moves carry the F-walk; muca-on-F in accept | 4–6 wk (≈95% new module; incidence/sign + snake estimator are the load-bearing risks) | Highest: unlocks the deferred matter line + ΔF(L) + σ₁/|P₁| via snake ratios; needs g(A,B) gauge muca + measured RT for full certificate | 8 / 7 / 7.5 | **7.5** | **PRIMARY.** Only design ≥7 on all lenses; failure modes absent by construction. Gate L=16/high-q on round-trip tests. |
| 2 | **Accept-side muca-in-B (Metropolis)** | YES for moderate κ / L≤8; transport (D_B) collapse + extensive-V cost are the risk at deep κ | 1.5–2.5 wk (reuses MucaB + U1LLR site block; sign-convention + invariant-distribution care) | Delivers canonical P(B), equal-weight κ_c, ΔB; right tool for the *moderate*-κ triple-point neighborhood | 6.5 / 7 / 6.5 | **6.7** | **FALLBACK + INSURANCE.** Cheap, share-everything with dual analysis; build in parallel. Bias EVERY B-changing kernel; freeze g for production. |
| 3 | **Bias-ladder (multicanonical PT)** | Conditional: cost scales with barrier L^(D-1) not latent-heat·V, but transport/σ_I unmeasured; 1D-g(B) non-ergodic at the 2D triple point | 8–11 wk (effort understated; U1HMC/U1LLR struct mismatch; thread-over-replicas fights inner OMP) | Same matter-line payoff as #2 but pricier; M~ΔF/2 could blow past 32 cores at q≥5 | 6.5 / 7 / 6.5 | **6.7** | **DOWN-RANK.** Strictly more complex than #2 for the same deliverable; revisit only if #2's V-scaling wall bites and the dual stalls. |
| 4 | **A-window LLR + matter overrelax + replica-exchange** | **NO** — crosses the (already-solved) A-jump and *assumes* B is slaved to A; adds no non-local move in B; stitched P(B) biased at coexistence cell | 7 (spine real, but exchange manager from-scratch; radial heatbath from-scratch; "embarrassingly parallel" false) | Line-tracer only; order ΔF_B untrustworthy; no phase labels | **4** / 7 / 6 | **5.7** | **KILL as a tunneler.** Right surgery (drop the B-box) on the wrong axis. Salvage only the radial-overrelax kernel as a decorrelator for #1/#2. |
| 5 | **Cluster / reflection / embedded updates** | **NO** (self-admitted) — reflections preserve \|φ\|; Gore–Jerrum torpidity; κ=∞ "anchor" physically empty (matter eaten) | 4.5 / 7 / 4 | Decorrelation force-multiplier only; C0 phase-overrelax is exact+free; the real κ=∞ anchor is a gauge driver needing no cluster | 4.5 / 7 / 4 | **5.2** | **KILL as a tunneler; HARVEST C0/C1.** Keep the exact B-invariant phase overrelaxation (C0) and radial reflection (C1) as cheap adjuvants inside #1/#2. |
| — | **Order-FSS + κ→∞ matching methodology** | N/A — analysis + routing layer | 6–8 d code | Required regardless: turns sampler output into the publishable certificate; pure-Z_q driver + dB=0 hop are load-bearing | 7 / 8 / 6.5 | **7.2** | **ADOPT WHOLESALE.** This is the spine of Sections 5–6; its tunneling claims correctly limited to "route around the crux." |

**Bottom line:** #1 (dual) is the primary weapon; #2 (accept-side muca) is cheap insurance and the
moderate-κ triple-point tool; the matching methodology + pure-Z_q driver + hop move are built
regardless; #3 is parked; #4 and #5 are killed as tunnelers but donate two exact decorrelation
kernels (C0/C1) and the "drop the hard B-box" lesson.

---

## 4. Recommended architecture

### 4.1 Primary sampler — dual flux representation (`U1Dual<D>`)

**The exact map (no Villain needed — our cosine action dualizes via Bessel weights):**
- Gauge: `e^{β cos θ_p} = Σ_{n_p∈ℤ} I_{n_p}(β) e^{i n_p θ_p}`.
- Matter: double-expand each link's hopping exponential → net flux `k_ℓ = j−j̄ ∈ ℤ`, dimer
  `l_ℓ = min(j,j̄) ≥ 0`.
- Link constraint from ∫dθ: **`q·k_ℓ + [dn]_ℓ = 0`** (the charge-q generalization is the single factor
  of q), so `k_ℓ = −[dn]_ℓ/q` is derived and configs require divisibility `q | [dn]_ℓ`.
- Site constraint `∇·k = 0` (auto-satisfied; kept as a runtime invariant).
- Radial site weight `W(f_x) = ∫_0^∞ dr r^{f+1} e^{−r²−λ(r²−1)²}`, `f_x = Σ_{ℓ∋x}(|k_ℓ|+2l_ℓ)`.
- `Z ∝ Σ_{n,l} [Π_p I_{n_p}(β)] [Π_ℓ κ^{|k_ℓ|+2l_ℓ}/((|k_ℓ|+l_ℓ)! l_ℓ!)] [Π_x W(f_x)]` — **all weights
  positive, no sign problem.**

**Reaction coordinate:** `F = Σ_ℓ (|k_ℓ|+2l_ℓ)`, with **⟨B⟩ = ⟨F⟩/κ exact** (the only κ-dependence in
the weight is `κ^{|k|+2l}` per link). The symmetric→Higgs F-walk is carried by the **free,
q-independent dimer (l) moves** on the k=0 vacuum — this is *why* it tunnels for all q.

**Moves (integer Metropolis, log-space accepts, no force/integrator/window):**
1. l-move `l_ℓ → l_ℓ ± 1` (carries the F-walk).
2. cube-move (shift 6 faces of a 3-cube ±1) — preserves `[dn]`=k on every link (∂∂=0).
3. plaquette±q matter-loop (`n_p→n_p±q`, `k→k∓1` around ∂p).
4. winding/plane moves (topological completeness; exponentially suppressed).
5. (Phase 2, contingent) surface worm (Mercado–Gattringer–Schmidt arXiv:1211.3436).

**Muca-on-F:** add `g(F)` to the **accept ratio** of integer moves (ΔF∈{0,±1,±2} → bias is a tiny
O(1) ratio of tabulated numbers). Reuse `MucaB` (`src/u1/u1_mucab.hpp`) **verbatim as the weight
container** — `gval`/`wl_record`/`flatness`/`halve_f`/`save`/`load` — **never `gprime`** (the froze
path). 1/t schedule (Belardinelli–Pereyra); **freeze g for production**, reweight
`⟨O⟩=⟨O e^{−g(F)}⟩_bias/⟨e^{−g(F)}⟩_bias`.

**New modules:**
- `src/u1/u1_dual.hpp` — `U1Dual<D>`: int `n[n_plaq]`, `l[vol·D]`, derived `k`, `f[vol]`
  (incrementally maintained); log-space tables `ln I_n(β)`, `ln W(f;λ)`, ln-factorials (lazy
  extension; **use `std::cyl_bessel_i` for `I_n`**, log-sum-exp guard only for overflow); the 4 local
  moves + `sweep()`; Bessel-ratio A estimator, F/F-histogram, W-ratio |φ|²; constraint checkers.
- `src/u1/u1_dual_incidence.hpp` — **build and test FIRST**: `Plaq<D>` (plaquette index
  `pidx(s,μ<ν)`, link→plaquette incidence list with signs **locked to `plaq_angle<D>` orientation,
  `u1.hpp:23`**, the `[dn]_ℓ` curl). This is the classic dual-code bug farm; the codebase has zero
  incidence scaffolding (`geometry.hpp` is site-only).
- `src/u1_dual.cpp` driver (auto-discovered by Makefile `DRIVERSRC` glob): modes
  `run`/`validate`/`mucabuild`/`mucarun`/`defect`.
- `src/test_u1_dual.cpp` + `test/test_dual_incidence.cpp` — validation ladder (§6).

**Defect-line σ₁/|P₁| (the Z_q blind-spot fix):** a charge-m external line shifts the link constraint
to `q·k + dn = −m` along the contour; `⟨W_m⟩ = Z_defect/Z` by **snake-ratio chains** (each step an
O(1) ratio → exponential error reduction). For q≥2 a charge-1 line cannot be screened
(1 ≢ 0 mod q). **De-risk on a 4⁴ tiny lattice against a direct partition-function ratio before
committing the full estimator** — it has zero code reuse and is the sole reason to prefer the dual
over #2 for the blind spot.

**Insertion mechanics:** purely additive, header-only + one driver; U1HMC and the 1746-test suite are
untouched; merge-safe.

### 4.2 Parallel-insurance sampler — accept-side muca-in-B (`U1MucaMetro<D>`)

Lift the single-site scalar Metropolis from `U1LLR::update_constrained1d`
(`src/u1/u1_llr.hpp:450`, site block ~507–522), **strip the box guard**, and change the accept to
the codebase convention (biased = `exp(−S + g(B))`):

```
accept iff  u < exp( +κ·dB − dS_pot + (g(B_run+dB) − g(B_run)) )
```

**Convention is load-bearing** — this matches `MucaB` header `exp(−S+g(B))` and the reweighter
`u1_mucab.cpp:74` `exp(−g(B))`; the *opposite* sign reinforces the barrier and diverges WL the wrong
way. Add a `g≡0` byte-identity test and a known-linear-`g(B)=cB` test (must shift κ_eff→κ+c) as the
first gate.

Critical correctness items the verifiers flagged:
- **Bias every B-changing kernel**, not just the matter Metropolis: gauge link moves change B, so the
  g-ratio must enter the gauge accept too (or keep the gauge HMC force at plain κ but ALWAYS apply
  the accept correction `dH −= g(B_f)−g(B_i)`, `u1.hpp:359`). A canonical-invariant gauge kernel
  composed with a muca-invariant matter kernel shares **no** stationary distribution.
- **Freeze g for production** (WL build breaks detailed balance); jackknife ΔF over independent
  production segments.
- Drive the matter sector *primarily* by the biased local moves; do not let an unbiased scalar HMC
  relax B back across the gap between sweeps.
- Eliminate B-drift: do not carry an incremental B_run through any dB=0 sweep; recompute B once per
  measurement.

Donate two exact, FD-checkable decorrelation kernels from the killed designs:
- **C0 phase overrelaxation** (`φ_x → e^{2i·arg(F_x)} conj(φ_x)`): preserves |φ|, B, A, S_pot exactly
  → microcanonical, accept≡1, legal in any window with no bookkeeping. Use full-resultant `arg(F_x)`.
- **C1 radial reflection** (`ρ' = 2ρ̂ − ρ` about the local quartic minimum): the only B-moving local
  matter move; **seed ρ̂ from a ρ-INDEPENDENT value (analytic cubic root, not Newton-from-current-ρ)**
  or the involution breaks and silently biases ⟨|φ|²⟩; include the `ρ'/ρ` measure Jacobian.

### 4.3 κ→∞ matching infrastructure (build regardless — near-zero crux risk)

- **`src/u1/zq_gauge.hpp` + `src/u1_zqgauge.cpp`** — standalone pure-Z_q gauge driver, int links
  `n∈{0..q-1}`, `θ=2πn/q`, Metropolis/heatbath + a WL/LLR-in-A mode (copy the `update_constrained1d`
  gauge half; discrete single-link flips change A by O(2D), windows hold trivially). All observables
  are angle-based (`wilson_loop`, `polyakov_abs`, `creutz_jack`, `monopole`, `photon_structure`) so
  feeding `θ=2πn/q` reuses them unchanged.
- **`src/u1/zq_hop.hpp`** — the exact **dB=0** link hop `θ_μ(x) → θ_μ(x) + 2πm/q` (m∈1..q−1),
  accept `exp(−β·dA)` from `link_staple` (`u1_llr.hpp:113`). `dB=0` is exact by algebra
  (`e^{iq·2πm/q}=e^{2πim}=1`), the abelian generalization of the GH_GUPD center flip (q=2 special
  case, which fully cured SU(2) κ=8). Wire behind env `U1_ZQHOP` in `u1_scan.cpp`, one hop sweep per
  trajectory. **Make the hop multi-link/correlated** to survive large β (q≥5 needs β up to ~5–6 where
  single-link `exp(−β·dA)` freezes). Honest scope: it restores **gauge-sector** ergodicity on
  large-κ slices; it does NOT and cannot cross the matter B-jump (it is exactly orthogonal to F/B).

### 4.4 Analysis pipeline (shared by all samplers)

- **`scripts/u1_order_fss.py`** (new ~250 lines) — per-L: Borgs–Kotecky equal-**weight** pseudo-coupling,
  Lee–Kosterlitz `ΔF(L)=ln(P_peak/P_valley)`, Challa–Landau–Binder double-Gaussian latent heat,
  Binder U4 minima, χ_max scaling. Ingest BOTH canonical ts (dual F-hist or muca P(B)) and LLR tiles.
  **Add a reweighted-histogram P(O) path to `reweight.py`** — the existing `u1_llr_maxwell.py` is
  hard-wired to LLR rect-cell reconstruction and `reweight.py` exposes only means/susceptibilities,
  not distributions; this is ~150–250 new lines, not a "20-line shim."
- **`scripts/u1_triple_point.py`** (new) — three pairwise line intersections + spread = systematic
  error box; CC slope `dκ/dβ` vs `−ΔA/ΔB`; optional 3-basin joint (A,B) histogram.
- **`scripts/zq_match.py`** (new) — `Δβ_c(κ,L)` with `c/κ` fit; σ₁(β)/|P₁|(β) curve overlays;
  q≥5 window convergence; q=1 + ρ_M-collapse negative controls.
- **Reuse:** `scripts/u1_llr_maxwell.py` (`two_peaks`/barrier core), `src/measure/reweight.hpp` +
  `scripts/reweight.py` (MBAR, conventions pinned to (A,−B)), `src/measure/autocorr.hpp` (τ_int→N_eff).

---

## 5. Order determination & κ→∞ matching plan

### 5.1 Phase fingerprints (Elitzur-safe; all wired)

| Phase | σ₁ (charge-1 Creutz) | \|P₁\| | ρ_M | m_γ (structure factor) | matter L_link/χ_link |
|-------|----------------------|--------|-----|-------------------------|----------------------|
| Confined | >0 | ~0 | high | massive | — |
| Coulomb | =0 | >0 | low | **massless** | — |
| Higgs (residual Z_q, q≥2) | >0 (Z_q tension) | ~0 | low (screened) | massive | ordered (χ_link peak) |

- m_γ ONLY via the gauge-invariant magnetic structure factor `R⁻¹=a+b·p̂²` (`photon_structure.hpp`);
  the `photon_mass.hpp` temporal correlator is INVALID and its `summary.csv` columns are inverted —
  do not plot them. m_γ needs L_s≥16 (geometric anisotropy 16³×L_t; p̂²_min: 0.586@L=8 → 0.15@16).
- Never use ⟨φ⟩ (Elitzur). For fundamental matter Higgs↔confinement are analytically connected;
  q≥2 (center-blind) gives a genuine transition and keeps σ₁ a valid confinement test at all κ.

### 5.2 FSS volumes & order discriminator

Order is an L→∞ statement; a single L cannot settle it. **First order ⇔ ΔF(L) GROWS** as
`2σ_I·L^(D−1)+c` over ≥3 volumes; flat/shrinking ⇒ crossover. Corroborate with latent-heat
persistence, Binder-minimum deepening, `χ_max ∝ V`.

- **Gauge line:** L=6,8,10(,12) via validated 1D-LLR-in-A. Feasible now. Pre-register
  "inconclusive below L=12" for weakly-first-order segments (compact-U(1) order took L~14–18).
- **Pure-Z_q anchor:** L=8,12,16 (discrete, cheap); q=3,4 strongly first order → safely settleable.
- **Matter/Higgs line:** L=8,12,16 via the **dual** F-histograms (primary) or L=4,6,8 via
  accept-side muca P(B) (moderate κ). **Gate L=16 on a measured round-trip test** (round-trips/sweep
  must not fall faster than ~L^−4 from L=8→12 extrapolation).

### 5.3 Run grid (planning only — needs lead go-ahead)

| Block | (β, κ, L) | Machine / cost |
|-------|-----------|----------------|
| Dual cartography | q=1,2,3,4,5,8; κ-sweep 0.3→8; β-grid across lines; L=8 (ms-scale sweeps) | lenore+lucia, days (embarrassingly parallel over (β,κ,q)) |
| Dual matter-line FSS | line points × L=8,12,(16 gated); muca-on-F | lenore; L=16 only at line-crossing points |
| Pure-Z_q canonical | q=2,3,4,5,6,8; β bracketing CJR transitions; L=8,12,16 | lucia, <2 days |
| Pure-Z_q WL/LLR-in-A anchor | q=2,3,4 at β_c; L=8,12,16 | lucia, ~1 day |
| Matching slices (U1+q, hop ON, cold start) | κ=2,4,8(,16 for q=2); β-scans; L=8 all q, L=12 q={2,5} | lenore, ~3–4 days |
| Triple-point region | q=2(→3,4): β∈[0.7,1.1]×κ∈[0.2,1.0] step 0.05, L=8; ridge L=12 | lenore, ~2 days; this is the moderate-κ regime where accept-side muca + shallow κ-PT (proven connected, 65/65 pairs accept) is viable |
| m_γ corroboration | 16³×8, ~6–10 pts straddling Coulomb–Higgs line, q=2,5 | lenore, ~1 day |

### 5.4 Exact paper-claim criteria

**"Finite-κ triple point exists for q≤4":**
1. Three locator lines, each traced with its own family, intersect pairwise within a common error box
   (jackknife + MBAR + L=8→12 shift ~1/V); box stable under volume.
2. Clausius–Clapeyron `dκ/dβ = −ΔA/ΔB` within 2σ on each first-order segment, with ΔA, ΔB from
   double-Gaussian/Maxwell peak separations (ΔB from the dual F-histogram or muca P(B) — this is the
   leg that *requires* a working tunneler).
3. Corroboration: 3-basin joint (A,B) histogram at the candidate (β_t,κ_t). If unavailable, claim
   from (1)+(2) and label (3) as not-obtained (the standard already used for the q=1 junction).

Order labels attached ONLY where ΔF(L) grows over ≥3 volumes.

**"κ→∞ extrapolates to the discrete Z_q diagram":** a measured convergence between two simulations at
identical (L, β-grid, observables):
- `Δβ_c(κ,L)` decreasing, consistent with `c/κ`, κ=8(or 16) within 2σ of pure-Z_q at the same L
  (both β_c from the SAME equal-weight LLR criterion so the estimator systematic cancels).
- σ₁(β), |P₁|(β) curve distances shrinking monotonically κ=2→4→8.
- For q≥5: the finite-κ Coulomb window converging to the pure-Z_q two-transition window ⇒ κ_t=∞
  (this IS the #26 demonstration); for q≤4 no window ⇒ Coulomb capped, finite triple point.
- Negative controls: q=1 gives σ₁≡0 (screened) at all κ; ρ_M collapses on the Higgs side while σ₁
  stays finite (residual is Z_q, not trivial).

---

## 6. Phased implementation plan (each milestone has a validation gate)

**Phase 0 — κ→∞ infrastructure (1 wk, no crux risk, start immediately, runs alongside Phase 1).**
- M0.1 Build `zq_gauge.hpp` + driver. **GATE:** Z_2 self-dual β_c = ½ln(1+√2) = 0.4407 (CI assert);
  CJR phase counts (q≤4 one first-order; q≥5 two transitions + massless window).
- M0.2 Build `zq_hop.hpp` (+ multi-link variant). **GATE:** standalone test asserts `|dB| < 1e-13`
  for q=2..8, m=1..q−1, random φ/θ; detailed-balance test; hop accept >few-% at the largest β in each
  matching slice (else sector-mixing claim fails where needed).
- M0.3 Empirical hop check: rerun the L=2 hot/cold metastability gap (q=4,5,8, κ=4) with `U1_ZQHOP=1`
  — the **gauge-sector** gap must close. **GATE:** this is the go/no-go for the matching slices.

**Phase 1 — dual primary, built incrementally behind gates (Phase 1a–1e, ~4–6 wk).**
- M1a Incidence layer first: `u1_dual_incidence.hpp` + `test_dual_incidence.cpp`. **GATE:** signed
  plaquette sum around every link reproduces the curl of a random integer n-field; `∇·k=0`; each move
  preserves `q|[dn]`. Make-or-break correctness contract — pass before any physics.
- M1b Weight tables + core `U1Dual<D>` + 4 local moves. **GATE:** invariants every N sweeps;
  incremental f/F vs recompute < 1e-9.
- M1c Validation ladder (mandatory before any physics):
  - **κ=0:** pure-gauge dual (cube moves) vs HMC ⟨plaq⟩(β) ACROSS β_c≈1.0 (literature ~1.01).
  - **β=0:** Z factorizes (k≡0) → semi-analytic l-sector + W(f) check, **including at large κ**
    (f~10³–10⁴, where the W(f) Laplace asymptotics must hold — this is the deep-Higgs regime of
    interest).
  - **HMC cross-check:** 4⁴ & 6⁴ at moderate (β,κ), q=1,2: ⟨A⟩,⟨B⟩,⟨|φ|²⟩ agree within errors.
- M1d Muca-on-F: hook `MucaB::gval` into the integer accept (1/t schedule; `gprime` compile-guarded
  OFF). **GATE:** flat F-histogram + ≥10 round-trips across the B-jump at a known bistable point
  (q=4, κ~4, L=8); **measure round-trip time vs V at L=4,6 → extrapolate L=8,16 (GO/NO-GO for L=16
  order-FSS and the surface worm).**
- M1e Defect σ₁/|P₁| snake: prototype on 4⁴ q=2 vs direct Z_defect/Z ratio FIRST. **GATE:** clean
  validation before the full driver; round-trip/defect-insertion diagnostic at (q≥5, small β, large
  κ) before claiming σ₁ is measurable there (the least-certain payoff).

**Phase 2 — accept-side muca insurance (1.5–2 wk, parallel to Phase 1, separate workstream).**
- M2.1 `U1MucaMetro<D>` (strip box guard, add g-ratio in accept, pin sign convention). **GATE:**
  `g≡0` byte-identical to canonical Metropolis; `g(B)=cB` shifts κ_eff→κ+c.
- M2.2 Bias every B-changing kernel; freeze-then-produce. **GATE:** reproduce the exact froze point
  (q=2, L=4, β=1, κ=0.3) with round-trips across the B-window; reweighted ⟨B⟩ matches direct HMC at
  a moderate-κ ergodic point.
- M2.3 Harvest C0 (phase overrelax) + C1 (radial reflection, ρ-independent ρ̂). **GATE:** C0 conserves
  A,B,S_pot to 1e-12; C1 involution test (apply twice → identity to 1e-12) + DB ratio test.

**Phase 3 — analysis pipeline (concurrent with Phases 1–2).**
- M3.1 `reweight.py` reweighted-P(O) path + `u1_order_fss.py`. **GATE:** validate ΔF(L)-growth
  analyzer on a known first-order control (4D Potts proxy / pure-Z_3,Z_4 — NOT Z_2, which is
  second-order self-dual and cannot validate a barrier-growth analyzer).
- M3.2 `u1_triple_point.py` + `zq_match.py`. **GATE:** reproduce the q=1 FS junction at
  ~(0.8485, 0.526) from the line-tracing + intersection machinery before any q≥2 headline runs.

**Phase 4 — production (only after the lead's go-ahead, samplers gated green).**
Order: pure-Z_q grids (lucia) → dual cartography (lenore+lucia) → matching slices (lenore) →
triple-point region (lenore) → gauge-line LLR-FSS → m_γ 16³×8 corroboration → q=3,4 triple-point
repeats → SU(2)→2T analogue pass on `gh_string.cpp` (GH_GUPD on; β_c(2T)=2.24(8) convergence +
rep-Polyakov screening fingerprint).

---

## 7. Risks, open questions, and de-scope order

**Top risks (with the gate that retires each):**
1. **Dual incidence/sign bugs** — the classic dual-code failure. Retired by M1a's incidence-test
   contract + runtime `q|[dn]`/`∇·k` asserts + the κ=0/β=0/HMC ladder.
2. **L=16 deep-Higgs round-trip cost** (equilibrium interface barrier is representation-independent).
   Retired by M1d's measured RT-vs-V gate; if it fails, L=16 order-FSS is de-scoped and the surface
   worm (Phase-2 dual contingency, +1 wk) is triggered.
3. **High-q σ₁ defect throttle** — the plaquette±q move acceptance ~`(β/2)^q/q!` is tiny at small β,
   and the snake at high q/small β is the headline blind-spot region. Retired (or scoped out) by
   M1e's dedicated round-trip diagnostic; decouple this from the bulk-B-jump tunneling claim (l-moves
   carry the bulk, this only throttles the *observable*).
4. **Accept-side muca invariant-distribution bug** (gauge moves change B) — retired by M2.2's
   "bias every B-changing kernel" + a control-point reweighting test.
5. **W(f) deep-κ asymptotics** bias exactly the regime of interest — retired by the β=0 large-κ
   semi-analytic check.
6. **Triple-point error box dominated by L=8→12 shift** for q=4 (CJR threshold, possibly large κ_t) —
   may force a weaker "κ_t(4) > X, Z_q-matching consistent" statement; ridge refinement at L=12(16).
7. **m_γ resolution** needs L_s≥16 anisotropy, separate from the scan grid; if skipped, the
   Coulomb–Higgs line rests on σ₁/|P₁| alone (defensible; a referee may ask for the photon).

**Open questions:**
- Does the dual's F-walk stay polynomial at L=16 / q=8, or does an interface/droplet barrier survive
  muca-on-F (→ surface worm mandatory)? *Decided by M1d.*
- Is the triple point reachable with the cheaper accept-side muca (moderate κ) so the dual's
  full-certificate effort can be reserved for the deep-κ matching only? *Decided by the M2.2 +
  triple-point-region pilot.*
- SU(N)→discrete H: the dual does NOT generalize (non-positive nonabelian character weights). The
  SU(2)→2T deep-κ crux remains open and still rides the GH_GUPD center-flip + a 2T-sector
  left-multiplication hop (to be built) + muca. Do not let the U(1) dual absorb the SU(N) timeline.

**De-scope order if time-limited (cut from the bottom):**
1. SU(2)→2T analogue pass (separate program; U(1) is the active campaign).
2. Surface worm (only if M1d shows local dual moves suffice at L≤12).
3. L=16 matter-line order-FSS (keep L=8,12; order then "evidence-grade," consistent with the locked
   "order deferred" scope).
4. m_γ 16³×8 corroboration (lean on σ₁/|P₁| for the Coulomb line).
5. q=3,6,7 dual cartography (keep q=1 control, q=2,4 for #28, q=5,8 for #26).

**Non-negotiable floor (what ships even in the worst case):** the **LOCATE** deliverable — three lines
traced with Elitzur-safe locators + the κ→∞ Z_q matching (pure-Z_q driver + dB=0 hop + σ₁/|P₁|
overlays) — is achievable *without ever crossing the deep-Higgs B-jump*, because the matching slices
stay inside the Higgs basin and the κ=∞ endpoint has no B coordinate. ORDER determination is the
stretch goal, gated on the dual's measured round-trip performance.
