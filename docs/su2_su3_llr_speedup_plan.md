# SU(2)/SU(3) LLR speedup plan (vetted, ranked)

Goal: cut the wall-clock of the LLR density-of-states scan that locates the first-order gauge-freezing
`beta_f(kappa)` for SU(2)->2T/2O/2I (and SU(3) later). An L=8 4-kappa scan is currently ~hours on a
32-core box; BO/BI are 2-4x costlier.

Every candidate below was scrutinized by two adversarial lenses (bottleneck-correctness and
novelty-realism). An item is ADOPTED only if it survived BOTH (real bottleneck, novel vs the prior
`docs/su2_hmc_speedup_plan.md` adoptions, not Amdahl-capped into irrelevance, not correctness-breaking).
Speedups quoted are the ADVERSARIALLY-REVISED numbers, not the proposer's headline.

## LLR cost anatomy (the Amdahl map everything is measured against)

Per kappa, two-sided annealing runs ~25 cell-solves (`su2_llr.cpp:121,124`). Each cell =
`nseed` seed trajectories + `rm_solve`'s `K*NRM` restrained trajectories (`llr_hmc.hpp:88,103-105`).
With driver defaults (nseed=40, K=40, NRM=40): per cell = 40 + 1600 = **1640 full GaugeHiggsHMC
trajectories, of which the RM `K*NRM` block is 97.6%**. So the RM trajectory COUNT is the whole
ballgame; every per-trajectory or per-RM-iteration win maps ~1:1 to wall.

Within ONE restrained trajectory the cost splits by rep:
- **SU(2)->2T/2O/2I (today's pain):** the matter ACCEPT-step `S_matter()` (= `pot.value()`, called
  twice/traj for Hi+Hf, `llr_hmc.hpp:74,77`) is ~50% of the trajectory and runs the UNCACHED
  `O(n_ch^2*n_gen*d^3)` Casimir path (`scalar_invariants.hpp:68`), while the MD force already uses the
  cached `dV_dphibar_cached` superop. The value path is the un-optimized twin of the already-optimized
  force path.
- **SU(3) (future):** dominated by the link exponential `expi<3>` (no closed form) and the gauge force;
  the matter side mostly stays on the cheap tensor path for the planned Sigma-group irreps.

---

## ADOPTED (ranked by leverage)

### 1. tau-adaptive / residual-gated K in `rm_solve`  ***[highest leverage]***
- **Change:** `src/hmc/llr_hmc.hpp:103-109` (the inner `for k in K` block, currently fixed K). Collect the
  A-series during `seed()` (40 already-equilibrated trajectories), estimate `tau_int` with the EXISTING
  Sokal estimator (`src/measure/autocorr.hpp:135`, `kSokalC=6`), and size per-cell
  `K = min(Kmax, max(8, ceil(c*tau)))`. Most cells are white-noise (`tau~0.5`) and need K~8-13; only the
  back-bending EDGE cells (`tau~6.7`) keep K~40.
- **Speedup (revised):** **~2.0-3.0x** standalone on the 97.6% RM fraction. Both lenses survived this;
  bottleneck-correctness landed ~2.2x (Sokal `c=6` floor -> mean K ~12.6), novelty-realism ~1.9x
  conservative / ~3.3x aggressive (`K>=4*tau` floor, mean K~11). Use **~2.2x** as the planning number.
  Do NOT credit the "2.8x combined" the proposer quoted -- that double-counted a separate NRM cut.
- **Effort:** moderate (the `tau_int()` machinery already exists; wire it to size K).
- **Correctness guardrail:** K is the BLOCK-AVERAGE/TRACKING window inside continuous RM (the one-time
  A->A0 relaxation is done by `seed()`, NOT the K-loop) -- so reducing K trades variance for cost and
  CANNOT bias the a1 fixed point (verified: RM-truncation bias is independent of K). Two real hazards,
  both mitigated: (a) floor `K>=max(8, 6*tau)` so high-tau edge cells stay well-sampled; (b) estimate
  tau from the seed window (post-transient), with a conservative multiplier, because a first-iteration
  tau is contaminated by the a1 transient. The existing `resid_e>0.04` filter in
  `scripts/su2_llr_betaf.py:28` is a backstop that drops any cell left under-sampled. Validate: the
  K-adaptive `a1(e)` curve and resulting `beta_f` must reproduce a fixed-K=40 reference on ONE kappa
  within stats before production. Touches no HMC detailed balance / MD reversibility / `|phi|=1`.

### 2. value-superop: cache `pot.value()` as a `d^2 x d^2` quadratic form
- **Change:** `src/action/scalar_invariants.hpp:68` (`CasimirChannels::value`) + `:206`
  (`MultiInvariantPotential::value`). Mirror the EXISTING force hoist `build_combined_superop`
  (`:102`): build `Q = sum_c f_c * Psuper_c` once in the ctor through the same `apply_proj` path, then
  `value()` becomes one quadratic form `vec(M)^dag Q vec(M)`, `M=phi phi^dag` (exploit rank-1:
  `V_c = phi^dag P_c(M) phi`). Add a `value_cached()` storing `Qsuper` alongside `Asuper`. Consumed via
  `scalar_action` -> `S_matter()` (`llr_hmc.hpp:42`, twice/traj).
- **Speedup (revised):** **~2.0x at BT/spin-3 (d=7)** (measured: S_matter 50.7% of a real trajectory;
  value() vs cached = 72-109x; projected traj 2.03x OMP=1 / 1.89x OMP=4), **~2.8x at BO/spin-4 (d=9)**
  (S_matter rises to ~65%), and larger again at BI (d=13). This is the un-optimized twin of the already-
  adopted force superop, so it is MULTIPLICATIVE with OMP and closed-form-D (disjoint slices). The
  proposer's 2.0/2.2/2.4x slightly UNDERSTATES the d-growth; real BO is ~2.8x. Rep-agnostic mechanism,
  so it also benefits SU(3) (upside).
- **Effort:** moderate (near-mechanical clone of `build_combined_superop` + a guard test).
- **Correctness guardrail:** `value()` enters ONLY the accept/reject H (same kernel at Hi and Hf), so MD
  reversibility/symplecticity, detailed balance, and the frozen `|phi|=1` sphere are untouched; `a1(e)`
  is RM-driven by the gauge action A (untouched) -> `beta_f` unbiased. The cached form REORDERS FP ops
  (not bit-identical), matching to projector precision (measured max rel err 6e-14 at d=7, 3.9e-10 at
  d=13 -- tighter than the 2e-8 the force superop already ships). Gate with an `fnorm` cross-check vs the
  uncached `value()` (same `enable_fast_if_valid` pattern), and confirm acceptance/`<A>` unchanged.

### 3. Two-phase RM schedule (find + average) in `rm_solve`
- **Change:** `src/hmc/llr_hmc.hpp:107` (the `a1 += 12/(W^2(m+1))*(meanA-A0)` step). The current harmonic
  `1/(m+1)` gain scales as `1/nplaq^2` and pre-decays before it can close a far bracket-to-cell gap: at
  L=8 it genuinely under-converges from the global `a1_cold=3.0`/`a1_hot=1.8` bracket. Split into a short
  FIND phase (slow-decay `g0/sqrt(m+1)`, L-ADAPTIVE gain -- not a fixed `g0=2x`, which diverges at L=4
  for high intensive susceptibility) then an AVERAGE phase with a fresh `1/(m'+1)` clock + trailing
  residual. **Also a one-line FREE rider:** `rm_solve` already accumulates the tail mean `sumA`
  (`:108`) but returns the raw last iterate (`:113`) -- return the tail-averaged a1 for noise reduction
  at zero cost.
- **Speedup (revised):** **~2.5-3.5x on the RM term** per the bottleneck lens (NRM 40->~20 is the solid
  ~2x; the larger figure needs the L-adaptive gain). The novelty lens was more conservative (~1.03x if
  the bias only bites the ~2 pass-start cells, since interior cells warm-start from the neighbor's
  converged a1). HONEST plan number: treat this as primarily a **CORRECTNESS fix** (removes a real
  under-convergence on the first cell of each pass + the back-bending region that `resid_e<0.04` does NOT
  catch -- the restraint pins `<A>~A0` while a1 stays far off, biasing `lnrho=cumint(a1 dA)`), with a
  modest **~1.2-1.5x** wall benefit on top of item 1. Do NOT bank the 3-5x; the K cut (item 1) and the
  schedule fix attack overlapping RM budget and must not be multiplied naively.
- **Effort:** moderate.
- **Correctness guardrail:** RM consistency requires `g*s in (0,2)` -> the find gain MUST be L-adaptive
  (Kesten/slope-normalized), NOT a fixed multiplier (the proposer's `g0=2x` diverges at L=4, the very
  scale it benchmarked). The averaging phase keeps `1/m` unbiasedness. Validate `beta_f` against a
  full-harmonic reference on one kappa at production L=8/high-kappa (larger autocorrelation) before
  adopting. a1 only sets `beta_eff`; HMC/MD untouched.

### 4. Per-cell flush of each row (crash/timeout recovery + live monitoring)
- **Change:** `src/su2_llr.cpp:105-113` (`run_cell` buffers into `rows`; print only at `:126-129` after
  BOTH passes). `printf` the row + `fflush(stdout)` immediately after `rm_solve` returns; keep the final
  sorted block as a redundant footer. **The `fflush` is load-bearing** (redirected stdout block-buffers
  ~70 rows >> 21 cells, so a bare printf still loses everything on a crash). Emit the SAME divided format
  as line 129 so `resid_e` stays in the units the `resid_max<=0.04` filter expects.
- **Speedup:** **0x compute** -- but the local box OOM-crashes under heavy MC, and this recovers an
  18/25-complete run instead of restarting from zero, plus enables aborting mis-tuned scans early on live
  `acc_hmc`/`resid`. Bill as robustness/wall-time recovery, NOT a throughput multiplier.
- **Effort:** trivial.
- **Correctness guardrail:** pure output reordering; `su2_llr_betaf.py::load` re-keys on `round(e,5)`,
  de-dups overlap by smallest resid, and re-sorts -- streaming unsorted rows is provably harmless. The
  inter-cell `a1g` carry is preserved (flush fires after `a1g=a1`).

### 5. SU(3): closed-form Cayley-Hamilton `expi<3>` (build BEFORE any SU(3) run)
- **Change:** `src/group/sun.hpp:11-29` -- add a `template<> Cmat<3> expi<3>` specialization next to the
  existing `expi<2>` at `:34`. Morningstar-Peardon: `exp(iH) = f0 I + f1 H + f2 H^2` for traceless
  Hermitian 3x3 H, with `f0,f1,f2` from the two invariants via the CH recursion. Call sites:
  `gauge_higgs_hmc.hpp:142` (drift), `hmc.hpp:50`.
- **Speedup (revised):** kernel **~2.3x** measured (NOT the proposer's 14-20x: at d=3 the matmuls are
  ~free, the cost is the transcendentals/cubic-root angle solve). Whole-trajectory: **~1.25-1.29x for
  FUNDAMENTAL SU(3)**, but only **~1.05-1.09x for the large-rep 2T/2O/2I SU(3) campaigns** -- the LLR
  engine ALWAYS runs the matter force (no gauge-only path), so for the headline Sigma(108)=(2,2)=d27
  irrep the matter `apply_tensor` (dimT=729) dominates and `expi<3>` is ~1-7% of the trajectory. Bank it
  as table-stakes SU(3) infrastructure that stacks with OMP, NOT as a needle-mover.
- **Effort:** moderate.
- **Correctness guardrail:** handle the degenerate-eigenvalue / small-argument branch (the f-coefficient
  formulas have a 0/0 at coincident eigenvalues -- use the MP Taylor-stable branch). Gate hard with a
  new `test_su3exp.cpp` vs the generic `expi<N>` (<1e-12; a CH reference matched to ~8e-16 at drift scale
  and unitarity/|det-1| ~1e-15). The generic `expi<N>` stays as the validated reference fallback.

---

## REJECTED (do not re-litigate)

| Candidate | Why rejected |
|---|---|
| **kappa-pass stream-parallel** (split cold/hot into separate procs) | ALREADY-KNOWN + bandwidth-saturated. kappa is one-per-process today (`su2_llr.cpp:56`) and production already runs `for k in 4 5 6 8` as 4 procs x OMP=7. The only new bit (2x pass-split) only pays through the sublinear, memory-bandwidth-bound OMP curve, which 28 threads already largely saturate. Real gain ~1.3-1.5x at best, often ~1x; the 22x was vs a serial strawman. |
| **kappa-process fanout "now"** | ALREADY DEPLOYED. This IS the live production config (`su2_bt_B/llr_prod_lenore.sh`). Incremental gain ~1x; it also targets the SMALL axis (M=4 kappas) not the ~20 cells. |
| **OpenMP cell-loop / cell-parallel (drop the carry)** | CORRECTNESS-BREAKING + subsumed. Two carries exist: the a1g warm-start (safe to drop) AND the neighbor CONFIGURATION carry that `seed()` never resets (`llr_hmc.hpp:88`) -- the actual annealing physics. Replacing it with one fixed cold/hot start forces the bistable MIDDLE cells (which set the Maxwell `beta_f`) to nucleate from a far start -> wrong-phase nucleation that `resid_e` does NOT reliably catch (a metastable config still pins `<A>=A0` with small residual). Also competes for the SAME cores as the already-adopted intra-trajectory OMP (all hot loops are already `omp parallel for`), so net incremental win is ~1-2.5x only at small L/few kappa, ~1x when the kappa axis already fills the box, with a real `beta_f` bias channel. |
| **MPI over cells** | ENVIRONMENT-WRONG. There is no cluster -- lenore/loranne/local are three independent single-node boxes with no interconnect/shared FS. A single 4-kappa L=8 scan exposes ~84 independent (cell,kappa) units >> 32 cores, saturable by trivial process fan-out. MPI buys ~0 on one node, adds launch/serialization/seed discipline, and existing `src/mpi/` is domain-decomposition (halo exchange), not a task launcher. |
| **adaptive-overlap-tiling (overlap 2->1)** | CORRECTNESS-BREAKING. Production data (`su2_bt_B/llr_L8_k*.out`) shows overlap cells are NOT redundant: in the back-bend ONE pass systematically FAILS to confine (resid 0.10-0.24, plaq pinned at 1.0) and the OTHER is the only usable point. Cutting overlap deletes a real confined data point at the low-e edge where the equal-area integral sets `beta_f`. Capped at ~8% even if safe. (The non-uniform-de half is speculative + chicken-and-egg.) |
| **rm-slope-init (extrapolated a1 guess)** | ~1.0x AS-IS. `rm_solve` is a HARDCODED fixed-iteration loop (`for m<NRM`, no early-stop/break). A better start changes only accuracy-at-fixed-NRM, NOT the trajectory count -> zero wall change. Its claimed win is 100% contingent on an unbuilt residual early-stop. (Fold the idea into the schedule fix if ever.) |
| **rm-polyak-ruppert as a throughput win** | AMDAHL-OVERSTATED. The pathology only bites the ~2 pass-start cells + back-bend (interior cells warm-start from the neighbor); a global NRM cut isn't unique to PR. Real wall gain ~1.0-1.15x, not 1.6x. The genuine value (a correctness fix) is folded into ADOPTED item 3; the constant-step `c` needs per-(L,kappa,lw) tuning and overshot in the proposer's own test. |
| **rm-skip-deep-branches** | AMDAHL-CAPPED ~1.1-1.3x + couples cells. Seed cost (40 traj/cell) is uncuttable, and this driver is a serial annealing chain (deep ends are where each pass STARTS), so under-resolving them feeds biased config/a1_0 into the middle. Dominated by item 1; only worth it as a residual early-stop refinement on top. |
| **fuse-A-into-gauge-force** | AMDAHL-CAPPED ~1.05-1.08x (fund), <1.05x on the hard 2T/2O/2I/SU(3) reps that motivate the hunt -- A() is rep-independent, the matter force grows with d, so the A-fraction shrinks exactly where you need speed. Also infeasible as a free read-out (`beta_eff` needs the GLOBAL A BEFORE the force scales by `beta/N`; staple-centric force vs plaquette-centric A double-counts), so it needs a real force-kernel rewrite for a single-digit-% fund-only win. Skip unless a large fundamental campaign appears. |
| **nmd cannot drop** | NEGATIVE RESULT (recorded, not a change). The harmonic restraint smooths only the GAUGE force; the stiff matter sector (frozen-phi geodesic + scalar) sets the integrator error. Deep cold/ordered cells collapse to ~0-3% acceptance for nmd<=20, recovering only at nmd>=26. Halving nmd would lose the cold pass that locates `beta_f`. DO NOT attempt an nmd cut. |
| **skip-value-in-set-beta** | 0x (informational). Confirms `value()` is already only in the 2 accept-step calls, not the kick loop -- it bounds item 2's win, contributes nothing itself. |
| **SU(3) Gell-Mann force projection** | AMDAHL-CAPPED ~1.02-1.03x. The 8x `trProd` projection is only ~6-10% of the gauge force, which is itself a shrinking slice for large-d reps. The "same idiom in `hop_link_g`" claim is wrong (that's a DMat matvec). MEDIUM sign/index risk for ~2% wall. |
| **SU(3) fund_alg guard** | TRIGGER NEVER FIRES. The planned Sigma-group breaking irreps (Sigma(108)={4,2}, Sigma(216)={5,1}, Sigma(648)={6,3}, Sigma(1080)={6}) ALL keep `use_fast=0` and never call `fund_alg<3>`. The proposer confused the physics label "(3,3)" with the literal Young-row irrep `{3,3}` (the only thing that flips `use_fast=1`), which is not in the program. Zero benefit; option (a) "force tensor" would be a 10.8x REGRESSION on `{3,3}` (tensor is slower there). |
| **cx-limited-range / drop-dense-zero-guard / SoA-rep / fno-math-errno / reuse-traj-buffers** | SU(2) DEAD PATH or rounding-error. The SU(2)->2T/2O/2I reps route through `su2_closed_D` (`rep_general.hpp:107`), which never touches `dmat_expi` or `DMat::operator*` -> those flags give EXACTLY 0x on the current campaign. `cx-limited-range` is a real ~1.5-2x SU(3)-only `dmat_expi` kernel win (~1.05-1.12x whole-LLR), worth banking for SU(3) but must be SCOPED (it regresses some closed-form paths and breaks bit-reproducibility). `fno-math-errno` ~1% (the libmvec-SIMD mechanism does NOT fire here -- verified scalar `sincos` symbol remains). `reuse-traj-buffers` targets ~0.01% of a trajectory (measured copy = 0.033%, NOT the claimed 1.02-1.05x; allocs are on the serial master thread, no OMP contention). |

---

## SU(3) readiness (build before any SU(3) LLR run)

1. **`expi<3>` Cayley-Hamilton specialization** (ADOPTED item 5) -- the only genuinely needed new closed
   form; gated by `test_su3exp.cpp`. ~1.05-1.29x depending on rep.
2. **`-fcx-limited-range` SCOPED to the SU(3)/non-closed `dmat_expi` TU** (NOT global -- it regresses the
   SU(2) closed-form inner ~7-18%). ~1.5-2x on `dmat_expi` -> ~1.05-1.12x whole-LLR for SU(3). Trivial
   effort; re-baseline reproducibility gold (changes last ULPs). Only worth it once an SU(3) rep with a
   non-trivial `dmat_expi` fraction is actually being run.
3. **value-superop (item 2) is rep-agnostic** -> carries to SU(3) for free; no extra SU(3) work.
4. The two adopted SU(2) engine speedups (closed-form `D^(j)(U)`, SU(2) closed `expi`/`fund_alg`) are
   N=2-only and do NOT carry over -- SU(3) needs its own closed forms; only OMP threading is free.
5. The catastrophic `fund_alg<3>=5820ns` general-link-log is NOT on the path for any planned SU(3) irrep
   (all keep `use_fast=0`); defer any fix until an irrep deliberately flips `use_fast=1` (none does).

---

## ORDER OF OPERATIONS (fastest path to cutting LLR wall-time)

For the SU(2)->2T/2O/2I campaign that takes hours TODAY:

1. **Per-cell flush (item 4)** -- trivial, ship first. Zero compute gain but immediately protects every
   long run from OOM/timeout loss and gives live monitoring. Do this before any other change so the
   validation runs below are themselves crash-safe.
2. **value-superop (item 2)** -- moderate effort, ~2.0x at BT / ~2.8x at BO, MULTIPLICATIVE with OMP and
   closed-form-D, low correctness risk (mechanical clone of the existing force superop + fnorm gate).
   This is the highest-confidence pure-compute win and the biggest single lever on per-trajectory cost.
3. **tau-adaptive K (item 1)** -- moderate effort, ~2.2x on the 97.6% RM trajectory count. Reuses the
   existing `tau_int()`; gate with the `K=40` reference cross-check on one kappa. This is the highest-
   leverage trajectory-COUNT cut.
4. **Two-phase RM schedule (item 3)** -- moderate effort. Primarily a CORRECTNESS fix (removes the
   under-converged pass-start/back-bend cells that `resid_e` misses) with a ~1.2-1.5x wall rider. Adopt
   the free tail-averaged-a1 return immediately; the L-adaptive find gain needs the production-L
   `beta_f` cross-check.

**Combined expected speedup (SU(2)->2T, honest):** items 2 and 1 are on DISJOINT slices (per-trajectory
cost vs trajectory count) and stack multiplicatively: ~2.0x x ~2.2x ~= **~4.4x** at BT. At BO (d=9) the
value-superop slice is larger (~2.8x) -> **~6x**. Item 3 adds a modest ~1.2-1.5x and a real bias removal
but overlaps RM budget with item 1, so do NOT multiply it in fully -- budget the realistic stack at
**~4-6x at BT, ~5-7x at BO/BI**, on top of the already-adopted OMP (OMP=4-8) and closed-form-D. Item 4 is
wall-time insurance, not a multiplier.

For SU(3): build `expi<3>` (item 5) + scoped `-fcx-limited-range` first (~1.1-1.3x), then the rep-agnostic
value-superop (item 2) carries over for free.
