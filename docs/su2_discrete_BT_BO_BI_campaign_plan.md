# SU(2) → discrete-nonabelian H (2T, 2O, 2I): unified staged frozen-Higgs (β,κ) campaign

**Status:** vetted master plan. Supersedes the three per-group drafts. Every load-bearing number/API below is grounded against the live repo this session; the per-group drafts' factual errors (μ² transcription, "χ_link must be wired", "5.82 is wrong dimension/unsourced", "D=3 screening is the Objective-2 primary", "no SU(2) gauge-mass operator exists") are corrected here and the corrections are flagged inline.

**Goal.** Break SU(2) → a discrete *nonabelian* subgroup H (BT=2T |H|=24, BO=2O |H|=48, BI=2I |H|=120) via a SINGLE large Higgs irrep, freezing the OVERALL SCALE |φ|=1 as κ→∞ (gauge digitization). Map the (β,κ) phase diagram and verify the κ→∞ residual little group is H — **from SU(2)+Higgs in 4D ALONE** (the load-bearing success criterion, `phase-diagram-plan-FS`).

---

## PART I — SHARED METHODOLOGY

### I.1 Frozen-sphere setup (identical for all three groups)

- **Action.** `S = -β Σ ReTr U_μν − κ Σ Re[φ†D^(R)(U)φ′]` (standard Fradkin-Shenker hopping sign), Wilson β=4/g², via `GaugeHiggsHMC<4,2>`. Driver `src/hmc_higgs_multi.cpp` (point driver: plaq, L_φ, L_link, χ_link, τ_int) and `src/gh_string.cpp` (adds σ_fund / Creutz / V(R) / fundamental Polyakov). Both honor `GH_FROZEN`, `GH_COLD`, `GH_GUPD`.
- **Freeze scale ONLY.** `GH_FROZEN=1` (`src/hmc/gauge_higgs_hmc.hpp`: `frozen_phi`, `project_pi_tangent`, `geodesic_drift`, `normalize_phi`) constrains |φ|=1; the **shape sphere stays fully dynamical** (great-circle geodesic drift, machine-precision constraint, exact+reversible — FD-verified at κ=8 and κ=64). **Never** build RATTLE / multi-constraint freezing — that would freeze the shape and mis-select H. Caveat carried for all groups: scale-freezing sends m_H→∞ (kills radial spectroscopy); fine for the phase map, out of scope here.
- **Potential = `MultiInvariantPotential` over the R⊗R̄ Casimir channels**, with the **validated per-group locking couplings** held FIXED across the whole (β,κ) grid so the residual H is fixed by construction. Use the COMPLEX rep (`GeneralRep<2>`, `real=false`); never `:real` — the `real_grad/real_hessian` path assumes a manifestly-real basis the tensor basis lacks and miscounts Goldstones (`docs/locking_couplings.md`).
- **Never pass `auto`.** `auto` sets all f_c=1 → reduces to (φ†φ)² = bare quartic = abelian isotropy (U(1)), NOT a discrete nonabelian H (`src/hmc_higgs_multi.cpp`). The explicit comma list is mandatory.

**Per-group fixed inputs (verified live: `build/singlets`, `build/test_align`, `docs/locking_couplings.md`):**

| H | spin j | rep `{2j}` | d | mult | μ² | mass gap | lock |
|---|---|---|---|---|---|---|---|
| 2T | 3 | `{6}` | 7 | 1 | **0.113** | 0.453 | RIGID |
| 2O | 4 | `{8}` | 9 | 1 | **0.108** | 0.126 | RIGID |
| 2I | 6 | `{12}` | 13 | 1 | **0.065** | 0.213 | RIGID |

> **CORRECTION (2I draft error):** the 2I draft wrote `mu2 = 0.213`; that is the **mass gap**, not μ². Use **μ²=0.065** (`docs/locking_couplings.md`). Better: read μ² from `test_align` / use the radial-stationary `radial_mu2` for |φ|=1 (for a frozen run the bare μ² is largely irrelevant once the radial mode is constrained out, but a wrong value mis-tunes the Stage-0 smoke).

**f_c (channel-C2 order printed by `./build/hmc_higgs_multi <2j>`):**
- **2T** (C2 `0,2,6,12,20,30,42`): `0.1287,0.1548,0.1835,0.2399,0.0056,0.1745,0.1130`
- **2O** (C2 `…,56,72`): `0.0562,0.2045,0.2230,0.2206,0.0444,0.0740,0.0066,0.0340,0.1366`
- **2I** (C2 `…,90,110,132,156`): `0.1182,0.0819,0.0353,0.0274,0.1550,0.0828,0.0157,0.0702,0.1083,0.1143,0.0189,0.1124,0.0595`

All center-blind (N-ality 2j mod 2 = 0): the Z₂ center (−I ∈ 2T⊂2O⊂2I) survives → the N-ality-1 **fundamental** Wilson loop is unscreenable by the Higgs at all κ. This is the structural fact the whole gauge axis rests on.

> **Tree-level caveat (incorporate from critiques C2/2O-A2/2I-C11a):** `find_stable_couplings` proves only that the **classical** Hessian at the singlet VEV is PSD with the right zero-mode count at ONE candidate point; `docs/locking_couplings.md` states the global-minimum-vs-other-subgroups question is *not* checked. "H by construction" is therefore a **classical tree-level** statement. The quantum lattice could in principle realign (Coleman-Weinberg), and spin-6 is exactly where 2T has a mult-2 singlet (a θ-modulus) competing with the 2I lock. → an **intrinsic ensemble check** (I.4 / the VEV-overlap order parameter) is REQUIRED, not optional. We do not assume; we measure.

### I.2 Observable set — valid discriminants and the BANNED traps

**Two INDEPENDENT axes. Classify every point on BOTH; never read the Higgs off a gauge observable.**

| Axis | Observable | Driver / location | Why valid |
|---|---|---|---|
| **GAUGE** (β axis): confined↔deconfined/freezing | **χ_plaq = V_plaq·Var(plaq)** (primary line finder) | **GAP — wire one line** `plaq.susceptibility(V_plaq)` into `gh_string.cpp` (`plaq` Stats accumulated at :159/:183, never emitted) and `hmc_higgs_multi.cpp`. **V_plaq = vol·D·(D−1)/2** (the plaquette count), NOT n_bonds — using the wrong V rescales the FSS peak HEIGHT (location is fine). | Cheap, sharp first-order-freezing specific-heat peak; the GLL freezing transition is a bulk plaquette discontinuity. |
| **GAUGE** | σ_fund / Creutz χ(R,R) / V(R)-fit / fundamental \|P\| | `gh_string` + `creutz_jack.hpp` + `string_tension_report.hpp` | N-ality-1 fundamental unscreenable at all κ. **Resolution floor:** `Rmax=min(L/2,4)` (`gh_string.cpp:112`) → R≤4 even at L≥12; σ(β) is only quotable at L≥12 and carries a ~17% χ-vs-V finite-T/L systematic the jackknife does NOT cover. Use σ_fund as a confinement *scale* only where the residual confines; deep in the Higgs region label by χ_link + \|P_fund\| (see I.3). |
| **MATTER** (κ axis): Higgs↔symmetric | **χ_link = V·Var(L_link)** (primary, ALREADY WIRED) | `hmc_higgs_multi.cpp:131-132`, `gh_string.cpp:208-212` | L_link is fully color-contracted ⇒ **gauge-invariant**; its Higgs-blindness is Fradkin-Shenker/Osterwalder-Seiler analyticity, **NOT Elitzur**. The ONLY matter-line separator. Peak grows with V ⇒ transition; flat ⇒ crossover. |

> **CORRECTION (2I/2O drafts):** χ_link is **already emitted** by both drivers — do NOT list "wire χ_link" as a Stage-0 task. The ONLY genuine wiring gap is **χ_plaq**.

> **Normalization hazard (critique m1):** the two drivers use different V for χ_link (`hmc_higgs_multi` uses vol; `gh_string` uses n_bonds=vol·D). Location is unaffected; **peak-HEIGHT FSS must use one driver's series only** — designate `gh_string` (n_bonds) as the FSS series for both χ_link and χ_plaq, and document V explicitly in every figure.

**BANNED discriminants (state these in the paper so reviewers don't ask):**
- **Helicity modulus Υ** (`src/u1/helicity.hpp`) — TESTED+REJECTED; it is a confinement OP (Υ≠0 in *both* deconfined and Higgs), redundant with σ_fund, U(1)-only, not wired for SU(2). Not the Coulomb-like demarcator.
- **FMS naive vector mass** φ†T^aD(U)φ′ — gauge-VARIANT (Elitzur). Do not use as an order parameter.
- **Tr(ΦF)-type general-irrep vector interpolator** `O_μν = Σ_a n^a F^a_μν`, n^a=φ†T^a_Rφ (`gaugeboson_op.hpp`) — **CORRECTION (2I/2O drafts over-claimed "no such operator exists"):** it DOES exist and is gauge-invariant to 1e-15 (the validated `Tr(ΦF)` generalization). For a discrete-nonabelian residual it returns **all-massive / no transverse massless pole** — so it is a legitimate *negative control* (confirms no continuous residual) but is **NOT a phase discriminator**, and a small measured mass does NOT prove discreteness. May be computed as a control; must not gate anything.
- "All gauge bosons massive ⇒ discrete H" — NOT a clean discriminator (weak-coupling-deconfined has near-massless gluons over a wide window).

### I.3 Objective 1 — "Coulomb-like": definition, locator, honest existence statement

**Definition (identical physics for all three groups; OUTPUT not target).** 2T/2O/2I are discrete-nonabelian ⇒ the κ→∞ residual has **no unbroken continuous U(1)** ⇒ **NO true Coulomb (massless-photon) phase** exists (Fradkin-Shenker criterion, `phase-diagram-plan-FS`). What exists is the **weak-coupling-deconfined region**: SU(2)'s asymptotically-free scaling window seen through the subgroup, where the discrete-confinement scale Λ is merely pushed exponentially far in β, terminating at the residual-H first-order freezing line.

**Honest framing of Objective 1 (resolving critiques 2T-M2, 2O-B1, 2I-C3 — they all converge):** this region is bounded ONLY by the gauge confinement/freezing line (toward confined) and the χ_link matter line (toward Higgs). It has **no internal order parameter and no third boundary of its own** — it is identified by *exclusion* (σ_fund below the L-resolution floor AND χ_link off-peak). The Fradkin-Shenker default null hypothesis is therefore that **"Coulomb-like" and Higgs are analytically connected** (no separating transition) for center-blind matter. So:

- **Objective-1 deliverable is a NULL result + two boundary lines, not a located fourth phase.** We report: (i) the gauge freezing line, (ii) the χ_link matter line, (iii) the quantified σ_fund sensitivity floor at our largest L (so "deconfined-by-resolution" has a number), (iv) the negative-control Tr(ΦF) confirming no massless mode. We do NOT draw an internal "Coulomb boundary" and we explicitly test (via χ_link FSS, Part I.5) whether any boundary between "Coulomb-like" and Higgs is real or the FS-expected crossover.
- **The global-U(1)_Φ Goldstone** (the +1 zero mode from the complex rep) is a *global* spectator phase, NOT a residual gauge U(1); it is invisible to σ_fund/χ_plaq/χ_link and must not be mis-read as a "massless mode" demarcating anything (pre-empts the G4 panic).

### I.4 Objective 2 — verify κ→∞ residual little group IS H (the RESOLVED method)

This is where the drafts erred most and the critiques converged hardest. The resolution below is binding for all three groups.

**The success criterion is fixed by `phase-diagram-plan-FS` (load-bearing, 2026-06-07): reproduce the KNOWN 4D discrete-group physics from SU(2)+Higgs ALONE, in 4D. "NO standalone H code and NO D=3 — both REJECTED as crutches." The freezing transition is the TARGET to reproduce, not an artifact to dodge.**

Consequently the Objective-2 spine is **the 4D large-κ gauge sector**, in this priority order:

1. **(PRIMARY, 4D) Reproduce the residual-H first-order freezing transition.** Scan β THROUGH β_c at large κ; show the GAUGE sector (χ_plaq specific-heat peak / σ_fund / \|P_fund\|) exhibits a **first-order** transition with confinement below it, and that **β_c(κ) saturates to a κ-independent plateau** as κ=4→8(→16). The plateau VALUE is an OUTPUT we quote. This is the only route that (a) satisfies "4D alone", (b) is not D=3, (c) is not question-begging.

2. **(PRIMARY, intrinsic, falsifiable) VEV-stabilizer overlap order parameter (NEW — closes the circularity hole, critiques 2T-C1/C2, 2O-A2, 2I-C11a).** On the sampled ensemble at large κ, project the gauge-invariant scalar bilinear onto the precomputed H-singlet ray (`vacuum_alignment.hpp::singlet_vev` gives the ray a priori) and confirm the overlap → 1 as κ grows, AND measure the broken-generator orbit-rank dynamically (is it 3 = full discrete break, or does it collapse to a continuous residual?). **G-INTRINSIC must be allowed to FAIL** — "residual is larger than H / a competing stratum / a continuous residual" is a real physics outcome (it is exactly the spin-2→Q8 soft-lock failure mode that already happened), NOT "a bug." This is the genuine, measured, falsifiable Objective-2 proof.

3. **(SECONDARY, group-resolving) 4D rep-resolved Polyakov/Wilson fingerprint** (`polyakov_loop_rep`/`wilson_loop_rep`, `observables.hpp:70/86`) across a probe-irrep battery at large κ. The PATTERN of which probe reps acquire nonzero loop / screen is the H character fingerprint and is computable in 4D (no D=3 needed). **The discriminating contrast (from `build/singlets`, verified live):** the *lowest* probe spin containing an H-singlet is **j=3 for 2T, j=4 for 2O, j=6 for 2I**. So probing j∈{3,4,6} and reading the contrast separates the three groups:
   - j=3 screens for 2T only (2O,2I: 0) → distinguishes 2T from 2O/2I;
   - j=4 screens for 2T,2O but not 2I → distinguishes 2O from 2I;
   - j=6 screens for all three (mult 2,1,1).
   A single "screened iff has singlet" pattern is NOT sufficient (critiques 2T-C1, 2O-2a, 2I-C7 correct); the **contrast across j=3/4/6 plus a Z₂-floor control** (every singlet-free center-blind probe must stay area/charged) is what earns the group label.

4. **(SECONDARY, D=3, OPTIONAL CROSS-CHECK ONLY) the screening table.** D=3 is where the residual confines and the area/perimeter contrast is sharp — but the master plan REJECTS D=3 as the success route. So the D=3 screening table is **demoted to an optional secondary fingerprint, off the critical path**, and must be labeled as such (resolving critiques 2T-A1/2O-FLAW1/2I-D1: the drafts wrongly made D=3 the Objective-2 primary). It is also contaminated by Higgs string-breaking at large κ (`discrete_observables_program.md`) and needs smearing+GEVP that **do not exist** in the code (`measure/smearing.hpp` absent) → at bare planar loops with `Rmax=min(L/2,4)` the j=6 (d=13) probe signal is noise-dominated by R=3. If run at all, it is qualitative.

**Stance on standalone-pure-H-MATCH vs intrinsic verification (the assigned tension, RESOLVED):**
- `no-modezq-endpoint-runs` (2026-06-16) forbids running a standalone endpoint and **using its value to window/anchor** the grid (confirmation bias). It is about grid placement and about not citing an imported number as a target/gate.
- `phase-diagram-plan-FS` (2026-06-07) settled that the residual is H **by construction** (irrep + tuned couplings), so a standalone pure-H run is **OPTIONAL POLISH**, not a proof requirement.
- **DECISION:** verify intrinsically (legs 1+2+3 above). A standalone pure-H gauge run is a legitimate *independent bug-check* (different implementation, clean β_c) ONLY if it never windows/anchors the grid — and since H is fixed by construction it is **not on the critical path**. The **dimensionless-ratio H-LGT match** (`discrete_observables_program.md` item H, comparing ratios in a common confining window away from transitions, lattice spacings need not agree) is the *cleanest* independent confirmation and is NOT forbidden by no-modezq (it does not place the grid) — we list it as a fundable follow-on, not a gate.
- **The literature freezing numbers (β_f: 2T=2.24, 2O=3.26, 2I=5.82, GLL arXiv:2208.12309/2312.10285, same Wilson β=4/g² convention — verified in `discrete-betac-hmc-cannot-tunnel`) are used ONLY as a post-hoc footnote comparison, NEVER as a pass/fail gate and NEVER to window the β-grid.** This resolves the internal contradiction the critiques flagged (drafts gated G-D/G5 on hitting 2.24/3.26/5.82, which is endpoint-anchoring). Window β GENERICALLY (uniform grid spanning the expected transition); quote the extrapolated plateau as output; compare to literature afterward.

> **CORRECTION (all three drafts):** the freezing β_f values are in OUR convention and ARE sourced (GLL, `discrete_observables_program.md:119,190`; 2.24(8) re-confirmed same-action in `discrete-betac-hmc-cannot-tunnel`). The 2I-critique-1 claim "5.82 is wrong dimension / unsourced" is itself wrong. But the FS-plan target is the *first-order freezing transition reproduced from 4D*, with β_f used only as a post-hoc check — exactly as above.

### I.5 How a (β,κ) scan is driven, and the SAMPLING REALITY (binding for all groups)

- **Engine:** `GH_FROZEN=1`, explicit f_c, COMPLEX rep. Pre-tune nmd per (L,κ) for accept≈0.65 — **start nmd LOW (~8–12) and tune UP**, do NOT start at 24 (frozen geodesic drift is exact; the integrator error is dominated by the gauge/coupled sector). Avoid per-point autotune (calibration-cost gotcha) — fix nmd per (L,κ). Consider **multi-timescale** (gauge coarse / matter fine) at large κ where the κ-stiffened matter force otherwise forces nmd up steeply.
- **The first-order freezing line CANNOT be tunneled by local HMC** (`discrete-betac-hmc-cannot-tunnel`, decisive 4-agent test). The apparent "β_c climbs toward β_f as κ grows" is the **hot branch trapping harder (spinodal divergence), NOT the transition moving.** Do NOT read β_c from raw hot/cold hysteresis.
  - **GH_GUPD center-flip** (`GH_GUPD="ncenter,nhb,nmetro"`) is now WIRED into BOTH `hmc_higgs_multi` and `gh_string` (committed). For center-blind reps the flip U→−U has dS_H≡0 exactly (verified 1.9e-14) → accepts on dS_g alone; the updater self-aborts loudly if the rep is not center-blind (`gauge_link_updater.hpp:85`). **Efficacy (measured):** FULLY cures κ=8 (hot=cold); at κ=32 it *crosses* the barrier but does NOT fully settle → **deepest κ still needs multicanonical.**
  - **GH_GUPD heatbath-veto leg** (`nhb`, N=2 Kennedy-Pendleton) decorrelates the gauge sector but **vetoes ~100% of proposals at deep κ** (the extensive matter-energy wall) → useful only at κ≲1. Use `ncenter≥1`; do not rely on `nhb` deep in κ.
  - **SU(2) multicanonical / Wang-Landau / LLR does NOT exist** (only `src/u1/*` U(1) versions). It is the prescribed tool at the deepest κ (`freezing_vs_nonabelian_endgame`) and is **unbuilt infrastructure (~300–500 lines, ~weeks).**
  - **BINDING RESOLUTION (resolves 2T-M4, 2O-A/B, 2I-1/C10):** the campaign **does NOT quote a single tunneled β_c at deep κ.** It (a) caps the routine center-flip-cured scan at **κ≤8** where the cure is verified, (b) locates the freezing line by **hot+cold spinodal bracketing**, reporting a **β window with stated systematic** (the 1980s hysteresis method), (c) reports β_c(κ) plateau saturation as the κ→∞ signature, and (d) lists "build SU(2) multicanonical-in-plaquette" as an EXPLICIT, separately-budgeted follow-on deliverable required *only if* a quotable single β_c at κ>8 is wanted. The cold/ordered branch is ntherm-stable at all κ; the bracket is the equilibrium-honest statement.
- **Autocorrelation is budgeted, not assumed.** τ_int(plaq) and τ_int(L_link) explode near peaks / on the freezing line (the U(1) analogue blew up ~80×). The drivers print τ_int/N_eff. **Gate: N_eff ≥ 200 at every line/peak point**; therm ≥ 20·τ_int; inflate errors by √(2τ). Run length floats UP to meet N_eff (this only makes the budget larger — see costs).
- **Starts:** use cold/annealed starts + hot-branch as the upper spinodal; never quote single-start error bars on deep-κ plaquette (within-basin blocking is meaningless when trapped).
- **Machines:** route ALL MC to lenore (32-core/125GiB). The local box is a **shared 15GB desktop currently at load ~21** — analysis ONLY (it OOM-crashed before). **No intra-job threading for SU(2)** (`GH_PARSWEEP` is U(1)-only) → parallelize across points (OMP=1 × many); large-L points are latency-bound single serial jobs with no speedup.
- **CI BLOCKER (promote from footnote, critiques 2T-m2/2I-D1):** the FD-force regression tests cover only the GeneralRep+MultiInvariant **link/matter-staple** force gap — the scalar-force path IS covered (`test/test_invariants.cpp`), but `test_scalar.cpp` link-force FD runs only Fundamental/Adjoint. Add a GeneralRep(spin-j)+MultiInvariant case to `test_link_force_fd` (reference `_review_tmp/fd_prod_check.cpp`). A silent regression here moves β_c and could change which stratum is selected → **this is a hard pre-campaign (G0) gate, not a closing note.**

### I.6 Honest cost model (replaces all three drafts' optimistic totals)

Anchors (measured, 4⁴ OMP=1): **2T 1.56, 2O 5.91, 2I 56 s/traj.** Per-traj scales as `≈ anchor · (L/4)⁴ · (nmd/nmd₀)`. **The drafts' "200–700 core-h / 1–1.5 days" totals are wrong** (they applied 4⁴ costs to 12⁴/16⁴ lines and omitted nmd and τ_int inflation; critiques 2T-1/2O-D/2I-3,4 are correct). Honest per-stage costs are rebuilt below from L⁴ scaling × explicit traj counts × N_eff inflation; **expect each group to run multi-day to multi-week on lenore, dominated by the L≥12 boundary/FSS points and the deep-κ τ_int blow-up.** 2I is ~36× the 2T per-traj cost (the d=13 fast path helps the hot loop but the 4096-dim basis is still built and every link touch is a 13×13 expm; memory-bandwidth bound → many OMP=1 jobs, never oversubscribe).

---

## PART II — PER-GROUP STAGE TABLES

Ordering BT → BO → BI (cheapest → most expensive). All stages: `GH_FROZEN=1`, explicit f_c, COMPLEX, pre-tuned nmd, cold+hot bracket, N_eff≥200 gate, route to lenore. Costs are honest order-of-magnitude lenore core-hours (OMP=1), nmd-and-τ_int-inflated.

### II.0 Shared Stage 0 (pre-campaign, ALL groups, BLOCKING)

| Goal | Action | Gate |
|---|---|---|
| FD-force on production path | add GeneralRep(spin-j)+MultiInvariant link-force FD to `test_scalar.cpp`; agree ~1e-8 | **G0a: BLOCKER** — no production until green |
| Wire χ_plaq | one line `plaq.susceptibility(V_plaq)`, V_plaq=vol·D(D−1)/2, in both drivers | G0b |
| Re-validate lock | `test_align <2j>`: mult=1, invariance<1e-6, Hessian PSD, nzero=3+1, no negatives; AND confirm a clear potential gap from the chosen H-stratum to the nearest competing-stabilizer ray (esp. 2I vs the 2T mult-2 spin-6 stratum) | **G0c** — if fail, re-run `find_stable_couplings` (budgeted separately, hours–day) |
| Smoke + per-traj timing | run one frozen trajectory each group; **record s/traj on lenore** before quoting any stage budget | G0d |

### II.1 BT = 2T (spin-3, d=7) — anchor ~1.56 s/traj @4⁴

| Stage | Goal | Mode | L | (β,κ) window | sampler/moves | ~core-h |
|---|---|---|---|---|---|---|
| A pilot | accept/nmd tune, frozen sanity (\|φ\|=1 to 1e-12), shape-sphere slow-mode τ_int check | hmc_higgs_multi/gh_string | 8⁴ | β∈{0.5,1,1.5,2,2.5}, κ∈{0,0.5,1}; +hot/cold pair at (1.0, 8) | HMC, nmd start ~10↑ | ~20–40 |
| B map | gauge freezing line (χ_plaq/σ_fund) + matter line (χ_link) across plane | gh_string + point | 8⁴ coarse, **12⁴** on lines | β∈[0.4,2.8] Δ0.2 (0.1 near peaks); κ∈[0,2.5] Δ0.25; **deep-κ β-line at κ∈{4,8} ONLY** (not 16) spanning β∈[0,3] | GH_GUPD `2,0,2` (center-flip) at κ≥1; hot+cold bracket | ~1.5–3k (dominant) |
| C intrinsic | VEV-overlap→1 + orbit-rank=3 at large κ; 4D rep-Polyakov j=3 (screens) vs j=2,1 (charged) | point + rep loops | 12⁴ | κ∈{4,8} at β below freezing | reuse B configs; cold | ~100–200 |
| D FSS/stitch | χ_plaq/χ_link peak height vs V (∝V ⇒ 1st order); β_c(κ) plateau; assemble 2-axis diagram | analysis | 8/12/(16 on lines) | on-line points only | √(2τ) inflation | ~200–400 (L=16 lines) |
| (opt) D=3 / H-LGT match | secondary screening fingerprint; dimensionless-ratio match | D=3 / hlgt.cpp | — | off critical path | — | follow-on |

**BT total ≈ 2–4k core-h (days–~1.5 weeks on 32 cores).**

### II.2 BO = 2O (spin-4, d=9) — anchor ~5.91 s/traj @4⁴ (~3.8× BT)

Same structure as BT; per-traj ~3.8× heavier, mass gap 0.126 (smaller → verify frozen-sphere acceptance/slowest-mode τ_int explicitly in Stage A, G0e). Stage-3 group-resolving contrast: **j=4 screens for 2O but j=3 does NOT (the 2O-vs-2T discriminator), j≤3 stay charged.**

| Stage | Goal | L | (β,κ) window | moves | ~core-h |
|---|---|---|---|---|---|
| A pilot | + frozen-sphere accept at the small (0.126) gap; heatbath-veto acceptance vs κ (expect useful only κ≲1) | 6⁴,8⁴ | β∈{0.8,1.2,1.6,2.0,2.6}, κ∈{0,0.5,1,2} | HMC, nmd↑ | ~60–120 |
| B map | freezing line + matter line; **deep-κ β-line κ∈{4,8} ONLY** | 8⁴, 12⁴ on lines | β∈[0.6,3.0] Δ0.1 GENERIC (do NOT center on 3.26); κ∈[0,3] Δ0.1 at β∈{0.8,1.2,1.6,2.0,2.6} | GH_GUPD `2,0,2`; hot+cold | ~3–6k (dominant) |
| C intrinsic | VEV-overlap→1, orbit-rank=3; 4D rep-Polyakov contrast j=3(charged)/j=4(screened) | 12⁴ | κ∈{4,8} | reuse B; cold | ~200–400 |
| D FSS/stitch | peak FSS, β_c(κ) plateau (compare to 3.26 post-hoc only), diagram | 8/12 (drop 16⁴) | on-line | √(2τ) | ~300–600 |
| (opt) D=3 / match | secondary | — | off path | — | follow-on |

**BO total ≈ 4–8k core-h (~1–2.5 weeks).**

### II.3 BI = 2I (spin-6, d=13) — anchor ~56 s/traj @4⁴ (~36× BT); the cost driver

Fast path ON (`exp(i w·T_R)`, self-check<1e-7) — but the 4096-dim basis is still built and each link touch is a 13×13 expm; memory-bandwidth bound → **many OMP=1 jobs, never oversubscribe.** Mass gap 0.213. **2T mult-2 spin-6 competing stratum** → G0c gap check is mandatory and the intrinsic VEV-overlap (Stage C) is the key falsifiable proof.

| Stage | Goal | L | (β,κ) window | moves | ~core-h |
|---|---|---|---|---|---|
| 0 | + re-validate 2I lock (was nsamples=8000); confirm 2I-vs-2T-mult-2 gap; FD on `{12}`; smoke timing | 4⁴/2⁴ | (β=2,κ=0.3) | — | hours–day |
| A pilot | accept/nmd; frozen sanity; NO κ≈0.3 *shape-sphere* slow mode (the real frozen failure mode, not a radial ridge) | 6⁴,8⁴ | β∈{1,2,4}, κ∈{0,0.2,0.5,1} | HMC, nmd↑, consider multi-timescale | ~0.5–1k |
| B matter line | χ_link(κ) peak vs L → Higgs line | 8⁴,12⁴ (FSS) | κ∈[0,3] at β∈{2,4} GENERIC; FSS L=8/12 | GH_GUPD `2,0,2` at κ≳1; hot+cold | ~3–6k |
| C gauge/freezing line | χ_plaq/σ_fund freezing line; **scan β THROUGH the transition at large κ, GENERIC window β∈[0,4] (do NOT center on 5.82)**; bracket hot/cold | 12⁴ (σ needs L≥12; honest R≤4 floor) | κ∈{0,1,2,4,8 cap} | GH_GUPD `2,0,2`; hot+cold spinodal bracket | ~4–8k (dominant) |
| D intrinsic | VEV-overlap→1, orbit-rank=3 at large κ (discriminates 2I from the 2T/2O strata); 4D rep-Polyakov contrast j=4(charged for 2I)/j=6(screened) | 12⁴ | κ∈{4,8} | reuse C; cold | ~0.5–1k |
| E stitch | β_c(κ) plateau (post-hoc vs 5.82 footnote), 2-axis diagram, σ_fund floor, Tr(ΦF) null control | analysis | — | — | analysis |
| (opt) D=3 / match | secondary fingerprint (NOTE: needs smearing+GEVP, unbuilt → qualitative only; j=6 d=13 probe noise-dominated at R≤4) | — | off path | — | follow-on (large-add) |

**BI total ≈ 8–15k core-h (multiple weeks); the 2I cost dominates the whole program.** Drop 16⁴ entirely; treat 12⁴ as the FSS ceiling, staggered early as long background serial jobs.

---

## PART III — GO/NO-GO GATES (per stage)

- **G0a (BLOCKER):** production FD-force on GeneralRep(spin-j)+MultiInvariant link path ~1e-8. NO-GO → fix before anything.
- **G0c:** `test_align` clean lock (mult=1, PSD, nzero=4) AND a clear potential gap to the nearest competing-stabilizer ray. NO-GO → re-search couplings.
- **GA (pilot):** accept 0.55–0.85 frozen; \|φ\|=1 to 1e-12; ⟨W[1][1]⟩==avg_plaquette to ~1e-15; **shape-sphere slowest-mode τ_int finite** (the actual frozen-HMC failure mode); heatbath-veto acceptance characterized vs κ. NO-GO → fix nmd/build/move-mix.
- **GB1 (gauge line):** χ_plaq & σ_fund features co-locate; **hot/cold bracket CLOSES at κ≤8** (cure verified there). At κ where the bracket does NOT close → that κ is beyond the center-flip cure; **report a hysteresis WINDOW, do NOT quote β_c**, and flag muca as the required follow-on. Never read β_c from raw hysteresis.
- **GB2 (matter line):** χ_link peak grows with V (transition) or flat (crossover) — **report whichever**; this is the arbiter of whether "Coulomb-like" and Higgs are separated (FS null = crossover/connected). N_eff≥200 required at the peak.
- **GC/GD-intrinsic (Objective 2, FALSIFIABLE):** VEV-overlap with the H-singlet ray → 1 as κ=4→8 AND dynamical orbit-rank = 3 (full discrete break). 4D rep-Polyakov contrast matches the group's j-ladder (2T: j=3 screens; 2O: j=4 screens, j=3 not; 2I: j=6 screens, j=4 not). **Allowed to FAIL** = "residual is larger/continuous/competing stratum" is a real result, not a bug.
- **G-stitch (Objective 2, 4D):** β_c(κ) saturates to a κ-independent plateau as κ=4→8; quote the plateau as OUTPUT. Compare to GLL β_f (2.24/3.26/5.82) ONLY as a post-hoc footnote; divergence-from-the-number is NOT a NO-GO (it is the expected finite-κ FS Higgs-transition behavior, not pure-H).
- **G-Coulomb (Objective 1):** the diagram shows exactly two boundary lines (freezing + χ_link) and the Tr(ΦF) control returns all-massive (no massless residual). If anything appears to be a "massless-photon Coulomb boundary" → artifact (no continuous residual exists) or the global-U(1)_Φ spectator; re-examine, do not publish a Coulomb phase.

---

## PART IV — KEY RISKS + MITIGATIONS (critiques folded in)

| Risk | Mitigation |
|---|---|
| **Deep-κ freezing barrier untunnelable by HMC** (the dominant risk) | Cap routine scan at κ≤8 (center-flip-cured); report hysteresis WINDOW + plateau, not a tunneled β_c; budget SU(2) muca-in-plaquette as an explicit follow-on only if a single κ>8 β_c is required. |
| **Objective-2 circularity** ("H by construction" is unfalsifiable; tree-level only) | Make the intrinsic VEV-overlap + orbit-rank order parameter the falsifiable primary; allow G-intrinsic to fail; verify the competing-stratum gap at G0c. |
| **Screening fingerprint doesn't separate 2T/2O/2I; D=3 contradicts the success criterion; needs unbuilt smearing** | Demote D=3 screening to optional secondary. Use the **4D rep-Polyakov j=3/4/6 contrast** (verified discriminating from `build/singlets`) + the 4D freezing-line reproduction as primary. |
| **Endpoint-anchoring** (gating on 2.24/3.26/5.82) | Window β GENERICALLY; quote extrapolated plateau as output; literature β_f is a post-hoc footnote, never a gate. |
| **Cost ~10× under-estimated; no SU(2) intra-job threading** | Honest L⁴×nmd×τ_int budget (Part I.6); drop 16⁴; stagger 12⁴ as background serial jobs; route to lenore, never the saturated local box. |
| **τ_int blow-up / hot-start non-equilibration at large κ** | N_eff≥200 gate; therm≥20τ_int; cold/annealed starts; √(2τ) error inflation; patch any copied driver (e.g. screening's hardcoded hot-start) to cold. |
| **σ_fund resolution floor (R≤4, L≥12, 17% χ-vs-V systematic)** | Quote the σ_fund sensitivity floor; do not present σ_fund≈0 as crisp; classify deep-Higgs by χ_link+\|P_fund\|, not σ_fund-as-scale. |
| **FD CI gap moves β_c / changes selected stratum** | G0a hard blocker. |

---

## PART V — PREDICTED PHASE STRUCTURE (OUTPUT, per group — falsify, do not target)

Common structure (two independent axes): a **gauge confinement/freezing line** in β (the residual-H **first-order** freezing wall, plaquette discontinuity, χ_plaq peak ∝V) and a **χ_link matter line** in κ. The deep-κ freezing β plateaus toward (post-hoc comparison) GLL β_f. The weak-coupling-deconfined ("Coulomb-like") region is the σ_fund-below-floor / χ_link-off-peak corner — **identified by exclusion, NO massless photon, NO third boundary**, and FS-expected to be analytically connected to Higgs (test via χ_link FSS). Regions seen on the two axes: confined-symmetric (small β, small κ); confined+condensed (small β, large κ — strongly-coupled-with-condensate; we drop the "residual-H confined" label there absent a group-resolved probe, critique m3); weak-coupling-deconfined (large β, small/moderate κ); Higgs (large β, large κ).

Per-group expected freezing scale (post-hoc check only, our convention): **2T ≈2.24, 2O ≈3.26, 2I ≈5.82.** Width-of-deconfined-band ordering across groups is NOT asserted (critique 2O-C2: the sign is not derivable a priori — the higher β_f opens the residual-confinement *later*, not necessarily a wider band at fixed β). Matter line likely a crossover at small β (the U(1)-frozen lesson: the link matter observable can be a center-blind crossover) — FSS decides.

---

## Key files (absolute)

- `/home/hlamm/Desktop/QC/higgs+gauge/src/hmc/gauge_higgs_hmc.hpp` (GH_FROZEN geodesic HMC)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/hmc/gauge_link_updater.hpp` (GH_GUPD: center-flip self-abort :85, heatbath-veto :150)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/hmc_higgs_multi.cpp` (point driver; χ_link :131-132; ADD χ_plaq)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/gh_string.cpp` (σ_fund/Creutz/V(R)/\|P\|/χ_link :208-212; ADD χ_plaq; Rmax=min(L/2,4) :112)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/action/vacuum_alignment.hpp` (singlet_vev, find_stable_couplings, radial_mu2)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/measure/observables.hpp` (wilson_loop_rep :70, polyakov_loop_rep :86, susceptibility :119)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/measure/gaugeboson_op.hpp` (Tr(ΦF) general-irrep vector interpolator — negative control only)
- `/home/hlamm/Desktop/QC/higgs+gauge/src/screening.cpp` (D=3 2T template; OPTIONAL secondary; 2T-hardcoded, spin≥5 OOM guard :78, hot-start :116 → must patch for 2O/2I + cold)
- `/home/hlamm/Desktop/QC/higgs+gauge/docs/locking_couplings.md` (f_c, μ²: 2T 0.113 / 2O 0.108 / 2I 0.065)
- `/home/hlamm/Desktop/QC/higgs+gauge/test/test_invariants.cpp` (scalar-force FD covered), `test/test_scalar.cpp` (link-force FD gap → add GeneralRep case)
