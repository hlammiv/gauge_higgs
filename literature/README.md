# Literature & Design Overview — Arbitrary-Dimension Nonabelian Gauge + Higgs HMC

This is the orientation document (a `literature/README.md`) for building an efficient Hybrid Monte Carlo (HMC) lattice code that simulates, in **arbitrary spacetime dimension D**, a **nonabelian gauge field** (Wilson plaquette action, group SU(N) — extensible to SO(N)/Sp(N)) coupled to a **Higgs/scalar field in an arbitrary irreducible representation R** (fundamental, adjoint, 2-index, …).

## 1. What you are building and why it is unusual

Almost every production lattice code is hard-wired to 4 dimensions and to SU(3) fundamental. The two hard requirements here — **arbitrary D** and **arbitrary irrep R** — are met by only two existing frameworks, and never both cleanly in one:

- **HILA** (CFT-HY) is the only surveyed framework with genuinely arbitrary *compile-time* dimension (`NDIM` macro) and a working gauge+complex-scalar Higgs app.
- **HiRep** (claudiopica) is the reference for arbitrary *representations*: its `autosun` build-time generator emits the represent/project maps for any (N, R), and — the key design insight — **the HMC force code never changes across representations because the force always lives in the Lie algebra su(N)**.

The synthesis target is: **HILA's compile-time `NDIM` + HiRep's representation-agnostic algebra-valued force + Grid's SIMD layout/expression-template lessons + a counter-based per-site RNG + openQCD's multi-timescale integrator algorithmics.**

## 2. The physics in one paragraph

The system has two field types living on different manifolds: gauge links `U_{x,μ} ∈ SU(N)` on links (a Lie-group manifold) and a scalar `φ_x` (a linear vector space in rep R, or a constrained sphere if the radial mode is frozen). The action is the Wilson plaquette action plus a gauge-covariant hopping term `−κ Σ 2Re[φ_x† D^(R)(U_{x,μ}) φ_{x+μ}]` plus an on-site potential `φ†φ + λ(φ†φ−1)²`. HMC augments both fields with Gaussian momenta, evolves a fictitious Hamiltonian with a symplectic integrator, and corrects all discretization error with a Metropolis accept/reject — making the algorithm **exact for any step size**. The single most important physics fact to reproduce is **Fradkin–Shenker**: for *fundamental* matter, the confinement and Higgs regions are **analytically connected** (no order parameter, only a first-order line ending at a critical endpoint or a crossover); for **trivial-N-ality matter (adjoint, etc.) a genuine phase transition can exist**.

## 3. Recommended reading order

**Tier 0 — Read first (the spine of the algorithm):**
1. Neal, *MCMC using Hamiltonian dynamics* (arXiv:1206.1901) — cleanest pedagogical HMC: detailed balance via reversibility + volume preservation, leapfrog, energy error, tuning.
2. Lüscher, *Computational Strategies in Lattice QCD* (arXiv:1002.4232) — the authoritative su(N) HMC derivation: algebra conventions, `H = ½(π,π)+S`, EOM `U̇ = πU`, the Wilson gauge force as TA-projection of link×staple, leapfrog, Metropolis. **This supplies most of the gauge-sector formulas verbatim.**
3. Duane–Kennedy–Pendleton–Roweth, *Hybrid Monte Carlo* (1987, no arXiv) — the original; short, worth reading once.

**Tier 1 — Gauge sector implementation:**
4. Morningstar–Peardon (arXiv:hep-lat/0311018) — exact SU(3) exponential via Cayley–Hamilton and the TA projection (reused for stout smearing later).
5. Takaishi–de Forcrand (arXiv:hep-lat/0505020) and Omelyan–Mryglod–Folk (arXiv:cond-mat/0110585) — the 2MN minimum-norm (Omelyan) integrator, λ≈0.1931833.
6. Sexton–Weingarten (1992, no arXiv) — nested multiple-timescale integration (split gauge vs Higgs vs fermion forces).

**Tier 2 — The scalar/Higgs sector and conventions:**
7. Csikor–Fodor–Hein–Jansen–Jaster–Montvay (arXiv:hep-lat/9507024) and Fodor et al. (arXiv:hep-lat/9409017) — the canonical SU(2)-Higgs action, the `Φ = ρ·α` decomposition, the bare-mass relation `a²m₀² = (1−2λ)/κ − 2d`, heatbath/overrelaxation/Bunk-reflection updates.
8. Bonati–Cossu–D'Elia–Di Giacomo (arXiv:0911.1721) — the **fixed-length (λ→∞) benchmark** and the warning that early "first-order" signals were finite-volume artifacts (Binder → 0.6667 crossover).
9. Günther–Hollwieser–Knechtli (arXiv:1908.10950) — concrete **constrained-HMC (RATTLE)** for gauge-Higgs; the reference when you actually run HMC on the scalar.

**Tier 3 — Arbitrary representations (the hard part):**
10. Del Debbio–Patella–Pica (arXiv:0805.2058) — building `U_R` and `T^a_R` from the fundamental: adjoint trace formula, two-index sym/antisym generators, Dynkin index/Casimir tables. The practical route to port.
11. Bertlmann–Krammer (arXiv:0806.1174) — generalized Gell-Mann basis for any SU(N), exact matrix elements.
12. Bonati–Franchi–Pelissetto–Vicari (arXiv:2106.15152) — cleanest modern adjoint lattice action `Tr[Φ^t U^adj Φ]`, gauge-invariant flavor order parameter.
13. Fonseca, *GroupMath* (arXiv:2011.01764) / Yamatsu (arXiv:1511.08771) — the highest-weight RepMatrices algorithm to port for truly arbitrary irreps, plus tables of dims/indices/Casimirs to validate against.

**Tier 4 — Phase structure & validation benchmarks:**
14. Fradkin–Shenker (1979, no arXiv) — the analytic-connection theorem; representation dependence.
15. Kajantie–Laine–Rummukainen–Shaposhnikov (arXiv:hep-lat/9510020) — the 3D electroweak benchmark, endpoint `x_c = 0.0983(15)`.
16. Hart–Philipsen–Stack–Teper (arXiv:hep-lat/9612021) and Kajantie et al. (arXiv:hep-lat/9811004) — SU(2)/SU(3) **adjoint** (Georgi–Glashow) phase diagrams; genuine transitions, `Tr Φ³` order parameter.

**Tier 5 — Software architecture:**
17. Boyle et al., *Grid* (arXiv:1512.03487) — virtual-node SIMD interleaving, expression templates.
18. Drach–Martins–Pica–Rago, *HiRep v2* (arXiv:2503.06721) + the HiRep repo — the representation-agnostic-force architecture.
19. Mazur et al., *SIMULATeQCD* (arXiv:2306.01098) — modern multi-GPU functor architecture, even-odd memory ordering.
20. Edwards–Joo, *Chroma/QDP++* (arXiv:hep-lat/0409003) — the layered data-parallel architecture lesson (and a cautionary tale on pre-C++11 template verbosity).

## 4. How the pieces fit

```
                    ┌─────────────────────────────────────────┐
                    │  HMC driver (Tier 0)                      │
                    │  heatbath p,π → MD trajectory → Metropolis│
                    └───────────────┬───────────────────────────┘
                                    │ calls abstract force()/exp()/drift()
        ┌───────────────────────────┼───────────────────────────┐
        ▼                           ▼                           ▼
 ┌─────────────┐          ┌──────────────────┐        ┌──────────────────┐
 │ Gauge sector │          │ Scalar sector     │        │ Integrator       │
 │ Wilson action│          │ hopping + V(φ)    │        │ leapfrog / 2MN   │
 │ staple,force │          │ rep R via D^(R)(U)│        │ nested timescales│
 │ (Tier 1)     │          │ (Tier 2/3)        │        │ (Tier 1)         │
 └──────┬───────┘          └─────────┬─────────┘        └──────────────────┘
        │                            │
        └────────► both forces project to su(N) (HiRep principle) ◄───┘
                                    │
              ┌─────────────────────┴──────────────────────┐
              ▼                                             ▼
   ┌────────────────────┐                       ┌────────────────────┐
   │ Group/Rep machinery │                       │ Lattice/Field infra │
   │ generators, exp,    │                       │ arbitrary-D geometry│
   │ D^(R), reunitarize  │                       │ even-odd, RNG, halo │
   │ (Tier 3)            │                       │ (Tier 5)            │
   └────────────────────┘                       └────────────────────┘
```

The decisive architectural lever is the **HiRep principle**: keep the integrator, momentum heatbath, and link update strictly in su(N); only the maps `R(U)` (group → rep) and `project` (rep d×d matrix → N²−1 algebra vector) differ between representations. Adding the Higgs requires *no new force machinery* — it adds a "matter staple" to the existing gauge-link force and projects with the same routine.

## 5. Pitfalls flagged repeatedly in the literature

- **Convention factors.** Papers differ on `½Tr` vs `Tr` in the plaquette, `κ` vs `γ` in the hopping, the `½` in the hopping, and Hermitian (`Tr T^aT^b=½δ`) vs anti-Hermitian generator conventions. Pick ONE, document it, and reproduce a source paper's *internal* consistency check before claiming a benchmark.
- **Elitzur's theorem.** Do not expect `⟨φ⟩≠0` as an order parameter for fundamental matter.
- **Finite-volume false first-order.** Early SU(2)-Higgs "first-order" signals at β≈2.5 were crossovers on large lattices. Always do volume scaling of the susceptibility peak and Binder cumulant.
- **Topological freezing.** Toward the continuum limit, autocorrelations of slow modes diverge (z≈5 for Q²); plan slow-mode-aware error analysis.
