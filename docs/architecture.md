# Recommended Code Architecture & Validation Plan

Greenfield C++ (C++17/20) project. Opinionated: combine **HILA's compile-time `NDIM`**, **HiRep's representation-agnostic algebra-valued force**, **Grid's SIMD layout + decltype expression templates**, a **counter-based per-site RNG**, and **openQCD's multi-timescale integrator algorithmics**.

## 1. Module / file layout

```
src/
  core/
    config.hpp            # compile-time D (NDIM), N (NCOL), Representation tag; build-time params
    geometry.hpp/.cpp     # arbitrary-D lattice: coord<->index, neighbor (hopping) table, even-odd masks
    field.hpp             # Field<T> over arbitrary-rank tensors; layout-abstract (AoS/AoSoA/SoA)
    parity.hpp            # EVEN/ODD/ALL subsets + X site iterator (HILA-style)
    rng.hpp/.cpp          # counter-based per-site RNG (Philox/sitmo); keyed on (seed, site, draw_idx)
    comm.hpp/.cpp         # D-dim Cartesian MPI; Cshift/Stencil halo exchange; overlap hooks
    expr.hpp              # ~small decltype expression-template engine (Grid lesson)

  group/
    sun.hpp/.cpp          # SU(N) matrix ops; exp (SU(2) closed, SU(3) Morningstar-Peardon, general scaled-Taylor)
    su2_quat.hpp          # quaternion SU(2) fast path
    algebra.hpp/.cpp      # su(N): T^a_F basis (generalized Gell-Mann), TA projection, project/reconstruct
    reunitarize.hpp/.cpp  # Gram-Schmidt + det fix; ‖U†U−I‖ diagnostic

  rep/
    representation.hpp    # abstract Representation: dim, T^a_R, N-ality, real/complex/pseudoreal, D_R(U), force-kernel
    rep_fundamental.hpp   # D_R(U)=U
    rep_adjoint.hpp       # D_adj(U)_{ab}=2Tr(T^a U T^b U†) (no exp); generators −i f^{abc}
    rep_twoindex.hpp      # 2S/2A generators + bilinear link
    rep_general.hpp       # highest-weight / weight-space engine (GroupMath RepMatrices port)
    rep_registry.hpp      # build-time (N,R) -> Representation; autosun-style generation/caching
    invariants.hpp        # T(R), C2(R), runtime self-consistency asserts

  action/
    gauge_wilson.hpp/.cpp # plaquette action, staple, gauge force F_g = (β/N)[U Σ]_TA
    scalar_higgs.hpp/.cpp # hopping (via D_R) + potential; scalar force; matter-staple link force F_H
    action_iface.hpp      # force(), value(), drift() abstractions consumed by the integrator

  hmc/
    momentum.hpp/.cpp     # Gaussian heatbath for P (su(N)) and π (scalar)
    integrator.hpp/.cpp   # leapfrog, Omelyan 2MN, nested Sexton-Weingarten, optional force-gradient
    constrained.hpp/.cpp  # RATTLE for frozen-length (|φ|=1) scalar
    hmc_driver.hpp/.cpp   # heatbath -> MD -> Metropolis; ΔH, accept; <exp(-ΔH)> & reversibility diagnostics

  measure/
    observables.hpp/.cpp  # plaquette, L_φ, L_link, Polyakov, χ, Binder, Tr Φ², Tr Φ³, Q^{fg}, monopoles
    correlators.hpp/.cpp  # FMS 0++/1-- correlators, smearing (APE/HYP), effective masses

  io/
    checkpoint.hpp/.cpp   # native binary: U, φ, P, π, RNG state, metadata (D, N, R, extents, couplings)
    ildg.hpp/.cpp         # ILDG/LIME read-write; converters to/from openQCD/MILC/NERSC

test/                     # unit + validation (see §6)
docs/                     # conventions.md (THE single source of truth for factors), theory notes
```

## 2. Data structures

- **Lattice geometry:** `D` is a compile-time constant. Site index ↔ D-tuple via mixed-radix; neighbor table `nbr[site][2D]` precomputed (forward+backward). All plaquette/staple loops are `foralldir(μ){ foralldir(ν>μ){...} }` over D — never hard-code 4. (HILA's `NDIM` model.)
- **Fields:** `Field<T>` templated on the per-site type (`SU_N` matrix for links, `RepVector<R>` or Hermitian matrix for the scalar, `Algebra` for momenta). Layout (AoS / AoSoA-AVX / SoA-GPU) is hidden behind `Field`; **never expose AoS site structs**. Use Grid-style virtual-node SIMD interleaving on CPU (over-decompose so boundary permutes are O(surface/volume)).
- **Gauge links:** `Field<SU_N> U[D]` (or one `Field` with a μ index). Stored once in the fundamental; rep-R action calls `D_R(U)` on the fly.
- **Momenta:** `Field<Algebra> P[D]` (store N²−1 real components for trivial Gaussian heatbath; convert to matrix when needed) and `Field<RepVector> pi`.
- **Representation object:** holds `dim`, `T^a_R` (cached), `N-ality`, real/complex flag, and two function pointers/virtuals — `D_R(U)` and `force_kernel` (the `iT^a_R` contraction) — plus `algebra_project`. This is the single swap point for arbitrary irreps.

## 3. Even-odd, RNG, communication

- **Even-odd:** store even sites before odd in memory (SIMULATeQCD). Provide `Parity` subsets (EVEN/ODD/ALL) for heatbath/overrelaxation sweeps and parallel decomposition. Not needed inside global HMC (all sites move at once) but essential for the local-update benchmark path and for halo bulk/boundary splitting.
- **RNG:** counter-based (Philox4×64 / sitmo) keyed on `(global_seed, encode(global_site), draw_index)`. Streams are a pure function of global coordinates ⇒ identical results regardless of MPI/GPU/SIMD decomposition; trivially parallel. Save/restore RNG state in checkpoints. (Prefer over Mersenne-per-rank or RANLUX-per-site.)
- **Communication:** D-dimensional Cartesian communicator, per-direction halos behind a `Cshift`/`Stencil` interface. Overlap: build the operator into the send buffer before sync, run halo exchange on a separate thread/stream, use GPU-aware MPI when present.

## 4. Backend abstraction

Functor + `RunFunctors`-style kernel iteration (SIMULATeQCD), selecting CPU-SIMD / CUDA / HIP / SYCL at compile time, with one centralized memory manager (named smart pointers). Physics code stays backend-independent.

## 5. HMC integration strategy

- Two momentum fields (su(N) link momenta + scalar momenta), both Gaussian-refreshed each trajectory.
- Total link force = `gauge_staple_force(U) + higgs_matter_force(φ,U)`, both ending in the same TA/`algebra_project`. **Unit-test each force separately by finite difference.**
- Nested integrator from day one: expose `n_gauge_steps`, `n_higgs_steps` (and later `n_fermion_steps`). Gauge on fine scale, Higgs on coarse. Default integrator = Omelyan 2MN (λ≈0.1931833); leapfrog kept for correctness checks. Force-gradient as an optional drop-in.
- Reunitarize links once per trajectory.
- Built-in diagnostics: `⟨exp(−ΔH)⟩ ≈ 1`; reversibility (forward then momentum-flip backward → return to start within ~1e-10 double); `⟨ΔH²⟩^½ ∝ ε²` (leapfrog/Omelyan) or `ε⁴` (FG) when ε halved.
- Frozen-length scalar (λ→∞): RATTLE constrained integrator (Lagrange-multiplier projection inside leapfrog) keeping `|φ|=1` to machine precision while preserving reversibility.

## 6. Staged validation plan (cheapest first)

**Stage 0 — Infrastructure unit tests.**
- Generator algebra: `[T^a,T^b]=if^{abc}T^c`, `Tr(T^aT^b)=½δ`; `f^{abc}` matches known SU(3) values.
- Rep invariants: `Tr(T^a_RT^b_R)=T(R)δ`, `T^aT^a=C₂(R)I`, dims/Casimirs vs Yamatsu tables; adjoint link orthogonal (`DD^T=I`, det+1).
- exp(iP) round-trips; reunitarize fixes injected noise.
- Force finite-difference test (gauge and scalar separately, all reps).
- RNG reproducibility across decompositions.

**Stage 1 — Z₂ gauge-Higgs (3D).** Cheapest physics validation. Ising link/site variables; check strong/weak-coupling and self-dual limits; target the multicritical XY/O(2) point of arXiv:2112.01824. Validates geometry, even-odd, observables, susceptibility/Binder machinery.

**Stage 2 — Pure SU(2) gauge.** Set κ=0; reproduce known average plaquette vs β (strong `P~β/4`, weak `P→1`). Validates Wilson action, staple, gauge force, su(N) exp, HMC accept/reject, `⟨exp(−ΔH)⟩=1`.

**Stage 3 — SU(2)-fundamental Higgs (frozen length, λ→∞).** Match pure-gauge plaquette at κ=0, turn on κ, reproduce the crossover/endpoint near **β≈2.72–2.73** (arXiv:0911.1721). Check Binder → 0.6667 (crossover) via volume scaling; confirm Fradkin–Shenker analytic connection (no order parameter). Cross-check the parameter map `μ₀²=(1−2λ)/κ−2D` against arXiv:2206.12884.

**Stage 4 — SU(2)-fundamental Higgs (variable length, finite λ).** Full HMC on the scalar (potential makes it nonlinear). Reproduce a published endpoint; verify integrator scaling and reversibility with both fields dynamical.

**Stage 5 — SU(2)-adjoint (Georgi–Glashow).** Switch the Representation object to adjoint (D_adj fast path). Reproduce the three-regime structure (SU(2)-confined / U(1)-Coulomb / U(1)-Higgs) of Hart–Philipsen–Stack–Teper (arXiv:hep-lat/9612021); monitor monopole density. Confirms the genuine transition expected for N-ality 0.

**Stage 6 — SU(3) and higher reps.** Generalize N: pure SU(3) plaquette vs β; SU(3)-fundamental Higgs; SU(3)+adjoint with the `Tr Φ³` order parameter (arXiv:hep-lat/9811004). Then exercise the general-irrep engine (2-index sym/antisym) with invariant self-checks. Optionally validate the 3D EWPT effective theory against `x_c=0.0983(15)` (arXiv:hep-lat/9510020).

**Cross-cutting at every stage:** measure plaquette, L_φ, L_link (primary transition locator), Polyakov, χ, Binder each config; do volume scaling to distinguish first-order (χ∝V, U4 dips) from crossover (χ saturates, U4→0.6667). Reproduce a source paper's *internal* consistency number before claiming a benchmark.

## 7. Build & IO

- CMake; compile-time `NDIM`, `NCOL`, `REPR`, `GAUGE_GROUP` (SU_N default; allow SO_N, quaternion-SU(2) fast path). An offline `autosun`-style generator (or `constexpr`/templates) emits per-(N,R) represent/project/exp routines so the HMC core never changes across reps.
- IO: ILDG/LIME read-write for gauge configs + native checkpoint storing U, φ, momenta, RNG state, and metadata (D, N, R, extents, couplings). Converters to/from openQCD/MILC/NERSC.
- Reference: HILA (`NDIM`, layout), HiRep (`autosun`, force-agnostic), Grid (SIMD/expr templates), SIMULATeQCD (GPU functors/memory/even-odd), Chroma (layering).
