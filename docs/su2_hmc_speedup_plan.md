# SU(2) Frozen-Higgs HMC Speedup Plan (BT / BO / BI campaign)

**Target:** `GaugeHiggsHMC<4,2>`, driver `build/hmc_higgs_multi`, frozen-sphere mode (`GH_FROZEN=1`).
**Reps gating the campaign:** BT = spin-3 ({6}, d=7), BO = spin-4 ({8}, d=9), BI = 2I = spin-6 ({12}, d=13).
**Status:** prioritized + adversarially vetted. Five candidates designed, prototyped in scratch dirs (production `src/` untouched), and measured end-to-end at matched physics. Headline at the bottom.

---

## 1. The measured bottleneck

Per-trajectory cost is ~100% **MD force evaluation**: linear in `nmd` (`cost/traj = a + b*nmd`, `b` = one Omelyan kick = one force eval); the per-traj fixed term `a` (2 `hamiltonian()` action evals + momentum refresh + save/restore) is small, and **freezing adds no cost**. The dominant force differs by rep:

| Tier | rep | d | fast path | dominant cost | breakdown |
|---|---|---|---|---|---|
| **BT** | spin-3 {6} | 7 | **OFF** (gated by `tensor_cost<=1.3*expm_cost`, rep_general.hpp:170) | `apply_tensor` D(U)φ (rep_general.hpp:196) = ~91-98% of kick | `scalar_force_hopping` rotate/rotate_dag ~59-63% + `add_matter_link_force` hop_link_g ~32% |
| **BO** | spin-4 {8} | 9 | **OFF** (same gate) | same `apply_tensor` D(U)φ | same split |
| **BI** | 2I spin-6 {12} | 13 | **ON** | `build_Dcache` -> `dmat_expi` (linalg.hpp:205, ~138 us/call, ~20 dense 13x13 matmuls/link, rebuilt every kick) = **84%** of kick | potential superop matvec #2 at 12-15% |

Key structural waste on the **tensor path (BT/BO)**: the same `D(U)·φ_y` is recomputed **uncached, 2x per link per kick** (`scalar_higgs.hpp:93` rotate vs `hop_link_g`'s `apply_tensor` at rep_general.hpp:154), and **~3x per link** across the three forces. The potential gradient `dV_dphibar` is already hoisted and cheap (5-7%); gauge force <2%.

**Measured anchors** (this box, OMP=1, frozen): L=8 spin-3 nmd=10 = **10.06 s/traj**. Clean s/traj L=6 nmd=8: fund 0.10 / spin-3 1.74 / spin-4 6.10 / spin-6 11.04.

Threading already exists on all three forces (`-fopenmp` is Makefile default); OMP=4 gives 2.0-2.7x under box contention, ~2.8-3.0x idle.

**Two integrator levers were tested and FAILED** (see Section 4): SW multi-timescale (matter force is both stiff AND expensive) and force-gradient 2MNFG (needs *more* steps at large kappa, 3 force-evals/step).

---

## 2. Prioritized candidate table

Speedups are **realistic, Amdahl-capped, end-to-end per-traj** factors measured at matched acceptance/`<exp(-dH)>` (not isolated per-kernel ratios). "Correctness gate" = the test that MUST pass before merge.

| # | Optimization | Realistic speedup (per-traj) | Impl cost | Correctness gate | Verdict |
|---|---|---|---|---|---|
| 3 | **OpenMP OMP=4-8** for gating large-L / 2I jobs (operational, zero code) | **2.8-3.0x @ OMP=4; ~5x @ OMP=8** (idle box). LATENCY only, not throughput | none | `test_scalar` 80/80 at OMP=4 AND OMP=8 (verified) | **ADOPT** |
| 1 | **SU(2) closed-form symmetric-power D^(j)(U)** drop-in (replaces tensor apply AND dmat_expi for all single-row SU(2) reps) | **BT 2.1x (nmd10) / 2.9x (nmd20); BO 3.7x; BI 1.35x** | medium | `test_scalar` 80/80 **AND** new `test_su2closed` 19/19 (cached==per-call BIT-IDENTICAL; reversibility spin-3/4/6; `<exp(-dH)>~1`) | **ADOPT-WITH-CARE** |
| 5 | **Eliminate duplicate D(U)φ_y** between scalar-hopping and hop_link_g (tensor path only) | **~1.4x BT/BO** (collapses ~3 D-applies -> ~2); no-op on BI | medium | `test_scalar` 80/80 (incl. force-FD spin-3/4/6) + BIT-IDENTICAL A/B for reps 6,8,12 | **ADOPT** |
| 2 | **Enable per-link D-cache for spin-4** (flip `link_cacheable` / lower fast-path gate) | **BO ~1.4-2.1x** (rises to ~2x at L=8/large nmd); **REGRESSION on BT** | low | force-FD + reversibility + `<exp(-dH)>` spin-4 {8}; assert spin-4 `use_fast==true`, spin-3 `use_fast==false` | **ADOPT (spin-4 only) / REJECT (spin-3)** |
| 4 | Force-gradient 2MNFG (OMF positive-substep 4th-order) integrator | **~0.5x (a SLOWDOWN)** in the large-kappa campaign regime | high | n/a (rejected) | **REJECT** |
| — | SW multi-timescale (prototyped earlier) | net-negative (0.268 vs 0.306 s/traj) | — | n/a | **REJECT** |

### Notes on the ADOPT-WITH-CARE / partial verdicts

- **#1 closed-form D^(j)(U)**: correct and exact (closed-form D vs tensor <1e-9 over 200 random links, unitary + homomorphism + cached==per-call bit-identical). The original "6-8x BT / 35x BO / 5-7x BI" claims were **D-build-in-isolation ratios** and are refuted at the traj level: real trajectories are **potential-bound** (`pot.value()` at d=13 is 86% of the BI Hamiltonian, ~77% for spin-3). End-to-end measured factors are the table values. Enabling it also turns the per-link D-cache **on** for spin-3/4, which is why it subsumes #2/#5 there.
- **#2 D-cache**: ADOPT for **spin-4 only**. The crossover (single `exp(iH)` build beats 3 tensor applies) is at d>=9. Spin-3 (d=7) is **below** it -> cache-ON is a **0.75x slowdown** (measured, refutes the "2-2.5x BT" claim). The gate must compare `3*tensor_cost` vs `expm_cost + 3*matvec` (NOT the current single-tensor heuristic, which would never flip spin-4 on).
- **#5 dedup**: subsumed by #1 (which enables the cache, eliminating the duplicate). It is the **fallback** if #1 is deferred — it is a smaller, lower-risk independent ~1.4x on the tensor path with proven bit-identity. Zero effect on BI.

---

## 3. Recommended implementation order

Biggest **safe** win first; lower-risk-first within a tier. Cumulative factors are multiplicative *only where independent* — they are NOT independent on the same rep (see caveats).

| Step | Action | Independent gain | Helps | Risk |
|---|---|---|---|---|
| **1** | **#3 OpenMP OMP=4-8** on the gating L=12/16 spin-3/4 jobs and ALL 2I/BI jobs. Set the env var; no build, no merge. | 2.8-3.0x (@4) / ~5x (@8) | **ALL tiers, esp. BI/2I** (scaling improves with volume: L6->L8 = 2.10x->2.94x) | none — pure operational |
| **2** | **#1 closed-form D^(j)(U)** drop-in. Port `scripts/rep_general_su2closed_PROTOTYPE.hpp` into `src/rep/rep_general.hpp`; add `scripts/test_su2closed.cpp` -> `test/`. | BT 2.1-2.9x, BO 3.7x, **BI 1.35x** | **ALL tiers**; biggest single code win for BT/BO. BI helped least (1.35x) — its bottleneck is the potential, not the rep. | medium — exactness gate must pass |
| **3** | **#2 D-cache flip (spin-4 ONLY)**. After #1, spin-4 already runs the fast path; lock the gate so spin-3 stays on tensor (or skip if #1 already enabled it). Mostly a *verification + gate-locking* step post-#1. | BO ~1.4-2.1x (already largely captured by #1's cache enable) | BO | low — but MUST lock spin-3 off |
| **4 (fallback)** | **#5 dedup** — ONLY if #1 is deferred. Independent ~1.4x on the BT/BO tensor path. Skip if #1 lands (the cache makes it a no-op). | ~1.4x BT/BO | BT/BO | medium |

### Cumulative expectation (realistic, honest)

- **OMP alone (Step 1)**: ~2.8-3.0x @ OMP=4, applied to every tier immediately, **no merge risk**. This is the single highest-leverage, lowest-risk lever and is operational *today*.
- **OMP @ OMP=4 x closed-form D (Steps 1+2)** — these ARE independent (one is parallelism, the other reduces per-thread FLOPs):
  - **BT (spin-3):** ~3x x ~2.1-2.9x ~= **6-9x** per-traj latency.
  - **BO (spin-4):** ~3x x ~3.7x ~= **~11x**.
  - **BI (2I spin-6):** ~3x x ~1.35x ~= **~4x** — the smallest gain, and BI is the most expensive tier (8-15k core-h). **BI is NOT unblocked by these rep-level changes.**

### Explicit note for BI / 2I (d=13)

The d=13 path is dominated by the **MultiInvariant potential** `pot.value()`/`dV_dphibar` (86% of the BI Hamiltonian), NOT the rep matrix. The closed-form D buys only ~1.35x; OMP buys ~3x. **The high-leverage next optimization for BI is the potential** (hoist/cache the channel superoperator in `pot.value()` and reuse the `dV_dphibar` structure — analogous to the unfrozen-2I hoisting in `gauge-higgs-hotpath-optimization`, which was NOT applied to the frozen Hamiltonian's `pot.value()`). That work is **out of scope of these five candidates** but is flagged as the real BI lever. Until it lands, **BI relies on OMP=8 (~5x) for latency** and raw core-hours for throughput.

---

## 4. Rejected ideas (do not re-try)

| Idea | Why rejected (measured) |
|---|---|
| **Force-gradient 2MNFG integrator** (#4) | **~0.5x slowdown** in the large-kappa regime where it was claimed to help. The frozen gauge-Higgs hopping force hits its per-step **stability** limit (large eps blows up regardless of order) before the 4th-order truncation term dominates, so FG needs *more* steps (nmd~20 vs ~13 for 2MN at acc~0.85), AND each FG step costs 3 force-evals vs 2. cost/accepted-traj: 2MN 2.77-2.81 vs FG 5.20-5.32 (spin-3); BI 5.41 vs 6.17-6.83. Implementation is correct (reversible to 1e-14, genuinely 4th-order at smooth points) — it just loses on this stiff constrained force. ADOPT only for a future smooth/small-step use case. |
| **SW multi-timescale (Sexton-Weingarten) integrator** | Net-negative: single-ts nmd=24 = 0.268 s/traj vs MTS nmd=6/ninner=2 = 0.306 s/traj at kappa=3. MTS pays off only when the stiff force is *cheap*; here the matter force is **both** the stiff force **and** the expensive force, so splitting it out gains nothing. |
| **D-cache for spin-3 (BT)** (the BT half of #2) | **0.75x slowdown.** At d=7 a single `exp(iH)` build costs *more* than the 3 redundant tensor applies it would remove (crossover is d>=9). BT keeps the tensor path. |
| **Helicity modulus** as a Coulomb/Higgs phase observable (separate program note) | Already TESTED+REJECTED — it is a confinement observable, redundant with sigma_1; not a speedup candidate but logged so it isn't revisited. |

---

## 5. Reusable benchmark harness (before/after verification)

All saved in `scripts/`, production `src/` untouched. Always measure **USER CPU** (`%U`), OMP=1, small jobs (L<=8) — the box is shared.

| Script | Purpose |
|---|---|
| `scripts/su2_hmc_bench.sh` | **Primary harness.** Builds clean + `-DGH_PROFILE` scratch drivers from unmodified `hmc_higgs_multi.cpp`; attributes per-traj cost via the `a + b*nmd` model. Subcommands: `build`, `nmd [L] [rep]`, `rep [L] [nmd]`, `vol [rep] [nmd]`, `prof [L] [rep] [nmd]` (component breakdown), `frozen [L] [rep] [nmd]`. |
| `scripts/gh_omp_bench.sh` | OMP scaling table (#3). Usage: `scripts/gh_omp_bench.sh 6 8 10 2 10 "1 2 4 8"`. |
| `scripts/gh_force_dedup_bench.sh` | Build + **BIT-IDENTICAL** A/B check + nmd-sweep timing for the dedup (#5). Usage: `scripts/gh_force_dedup_bench.sh src <prototype_src> 6`. |
| `scripts/gh_integrator_bench.sh` | cost/accepted-traj frontier for 2MN vs FG (#4 rejection evidence). |
| `scripts/bench_repapply_su2.cpp`, `scripts/bench_repD_su2.cpp`, `scripts/bench_kick_dcache.cpp` | per-kernel D-apply / D-build / kick-with-cache micro-benchmarks (#1/#2). |

### Prototypes + exactness gates (ready to port)
- `scripts/rep_general_su2closed_PROTOTYPE.hpp` — the #1 closed-form D drop-in for `src/rep/rep_general.hpp`.
- `scripts/test_su2closed.cpp` — the **new regression gate** (19/19) that MUST be added to `test/` (closes the gap that existing `test_scalar.cpp` runs reversibility/`<exp(-dH)>` ONLY on fund/adj, never on GeneralRep).
- `scripts/gh_fg_integrator_test.cpp` — FG correctness gate (kept for the record; FG is rejected).

### Standing regression gate (must pass before ANY merge)
`make build/test_scalar && ./build/test_scalar` -> **80/80**, plus (after #1) `./build/test_su2closed` -> **19/19**. Run at OMP=1, 4, AND 8.

---

## Headline

**Implement #3 (OpenMP, OMP=4-8) first — it is operational TODAY with zero code change and gives ~2.8-3.0x at OMP=4 (~5x at OMP=8) on every tier, with scaling that *improves* at the large volumes that gate the campaign.** Then merge #1 (closed-form SU(2) D^(j)(U)), gated on `test_su2closed` 19/19, for an additional ~2.1-2.9x (BT) / ~3.7x (BO) / ~1.35x (BI).

**Combined per-traj latency: ~6-9x (BT), ~11x (BO), ~4x (BI).** BI/2I (the 8-15k core-h tier) is the least helped because it is **potential-bound, not rep-bound** — its real next lever is hoisting/caching the MultiInvariant `pot.value()` superoperator at d=13 (out of scope of these candidates, flagged for follow-up). Reject the force-gradient and SW multi-timescale integrators (both measured slower on this stiff frozen force) and the D-cache for spin-3 (0.75x regression).

**Doc:** `/home/hlamm/Desktop/QC/higgs+gauge/docs/su2_hmc_speedup_plan.md`
