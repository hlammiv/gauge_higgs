# gauge+higgs — arbitrary-D SU(N) gauge + arbitrary-irrep Higgs HMC

A from-scratch C++ Hybrid Monte Carlo lattice code that simulates, in **arbitrary
spacetime dimension D**, a **nonabelian gauge field** (Wilson plaquette action,
group **SU(N)**) coupled to a **Higgs/scalar field in an arbitrary irreducible
representation R** of the gauge group.

Design follows the *HiRep principle*: the HMC integrator and force machinery are
representation-agnostic because the force always lives in `su(N)`; a representation
only supplies the rep generators `T^a_R` and the map `D^(R)(U)` from a fundamental
link to its rep matrix. D and N are compile-time constants (HILA/HiRep style);
the representation is selected at runtime.

See `docs/conventions.md` (single source of truth for factors/signs),
`docs/theory_notes.md` (implementation-grade derivations), `docs/architecture.md`
(full design + staged validation plan), and `literature/` (curated, arXiv-verified
papers + annotated bibliography).

## Build

Requires a C++20 compiler (g++ ≥ 11). No CMake/CUDA needed.

```bash
make                      # builds drivers + tests (NDIM=4, NCOL=2 by default)
make NDIM=4 NCOL=3        # rebuild for 4D SU(3)
make NDIM=3 NCOL=2 all    # 3D SU(2)
make test                 # build & run the full unit/validation test suite
make clean
```

`NDIM`/`NCOL` are compile-time (`-DNDIM= -DNCOL=`); the binaries are monomorphic and
fast. The test suite instantiates several (D,N) directly, so it is config-independent.

## Run

Pure gauge:
```bash
./build/gh_hmc  [L=8] [beta=2.3] [ntherm=50] [nmeas=100] [nmd=20] [tau=1.0] [seed=1]
# e.g. 4D SU(2) beta=2.3 -> plaquette ~0.603
```

Gauge + Higgs (arbitrary irrep):
```bash
./build/hmc_higgs <rep> <L> <beta> <kappa> <lambda> [ntherm] [nmeas] [nmd] [tau] [seed]
#   <rep> = fund | adj | <Young rows, e.g. 2,1>   (append ':real' for a real-scalar rep)
# e.g. ./build/hmc_higgs fund 6 2.3 0.20 0.5          # SU(N) fundamental
#      ./build/hmc_higgs 2:real 6 2.3 0.12 0.5        # SU(2) spin-1 (adjoint) via general engine
#      make NCOL=3 && ./build/hmc_higgs 2 4 5.7 0.1 0.5    # SU(3) sextet (d=6, N-ality 2)
#      make NCOL=3 && ./build/hmc_higgs 2,1:real 4 5.7 0.1 0.5  # SU(3) octet/adjoint (d=8)
# prints plaquette, L_phi=<phi^dag phi>, L_link (gauge-invariant hopping energy,
# the transition locator), acceptance, <exp(-dH)>.
```

Note: `NDIM`/`NCOL` are compile-time, so switching them rebuilds the drivers
automatically (a `build/.config-*` stamp tracks the active config).

MPI (domain-decomposed; `make mpi` builds both):
```bash
make mpi
# pure gauge:
OMP_NUM_THREADS=1 mpirun -np 4 ./build/gh_hmc_mpi 8 2,2,1,1 2.3
#   <L> <procgrid e.g. 2,2,1,1> <beta> [ntherm] [nmeas] [nmd] [tau] [seed]
# gauge + Higgs (any irrep):
OMP_NUM_THREADS=1 mpirun -np 4 ./build/gh_higgs_mpi fund 8 2,2,1,1 2.3 0.28 0.5
#   <rep> <L> <procgrid> <beta> <kappa> <lambda> [ntherm] [nmeas] [nmd] [tau] [seed]
#   prod(procgrid) must equal -np and each factor must divide L.
# Correctness: observables are identical for any procgrid (global-index RNG) — decomposition test.
```

## Source layout

```
src/core/    config, linalg (Cmat<N>/DMat), geometry (arbitrary-D), rng (counter-based),
             fields (gauge link + su(N) momentum), scalar_field (per-site complex d-vector)
src/group/   algebra (generalized Gell-Mann generators, TA projection, structure consts),
             sun (exp(iH): SU(2) closed form + general scaling-squaring, reunitarize)
src/rep/     representation (rep-agnostic interface + invariant checks), rep_fundamental,
             rep_adjoint (orthogonal link, no exp), rep_general (ARBITRARY irrep from Young
             rows / Dynkin labels via tensor-power + Young symmetrizer)
src/action/  gauge_wilson (plaquette, staple, force), scalar_higgs (covariant action, scalar
             force, matter back-reaction on links, observables)
src/hmc/     hmc (pure-gauge integrators + driver), gauge_higgs_hmc (combined)
test/        test_math, test_gauge, test_scalar   (1746 checks)
docs/        conventions, theory_notes, architecture
literature/  25 arXiv-verified PDFs + README (orientation) + BIBLIOGRAPHY
```

## Validation status (see docs/architecture.md §6 for the staged plan)

Every force is finite-difference-tested against the action; integrators are checked
for reversibility (~1e-9), 2nd-order energy-conservation scaling, and `<exp(-dH)>≈1`.

- **Stage 0 (infrastructure):** ✅ generator algebra, exp unitarity/det, rep invariants
  `Tr(T^aT^b)=T(R)δ`, `ΣT^aT^a=C2 I`, adjoint link orthogonality + homomorphism,
  gauge/scalar/matter force finite differences (SU(2), SU(3); fundamental + adjoint).
- **Stage 2 (pure SU(N) gauge):** ✅ cold plaquette=1; 4D SU(2) β=2.3 → plaquette 0.603
  (between strong/weak limits); reversibility, acceptance, `<exp(-dH)>`.
- **Stage 3 (SU(2)-fundamental Higgs):** ✅ finite-size validated — at β=2.3, λ=0.5 the
  confinement↔Higgs transition sits at κ_c≈0.265, and the susceptibility χ_link=V·Var(L_link)
  peaks there with the peak **growing with volume** (1.78 at 6⁴ → 2.42 at 8⁴) — the correct
  finite-size signature. Reproduce via `validate_su2higgs.sh` + `analyze_su2higgs.py`
  (→ `validation_su2higgs.png`). Reproducing a *published* endpoint number is the next step.
- **Stage 5 (adjoint / Georgi–Glashow):** ✅ correctness validated (forces, reversibility,
  `<exp(-dH)>`); quantitative phase diagram pending.
- **Stage 6 (general irreps):** ✅ `rep_general` validated — `[T^a_R,T^b_R]=if T^c_R`,
  dims/Casimirs, character cross-checks (general[1]=fund, [2]=adj, SU(3)[1,1]=conj-fund);
  SU(3) sextet/octet runs; adjoint fast-path vs tensor-engine agree at the simulation level.

## Roadmap

- **Performance:** OpenMP over sites ✅ done (~3-5× on 8-16 threads, volume-dependent);
  `GeneralRep` per-vector apply ✅ done (no d×d matrix; ~25× for large irreps). Next:
  thread-local scratch buffers for `apply_tensor` (kills per-apply allocation, the dim-48
  bottleneck); matrix-form (TA-projection) gauge force; specialize fundamental allocation.
- **Scaling:** MPI Cartesian decomposition + halo exchange — ✅ **done for the full
  gauge+Higgs stack** (`src/mpi/`, `make mpi`): `gh_hmc_mpi` (pure gauge) and `gh_higgs_mpi`
  (gauge+Higgs, any irrep). Validated bit-identical across np 1/2/4 (fundamental, adjoint,
  general reps) via the global-index RNG. Next: comm/compute overlap + fewer redundant halo
  exchanges; nested multi-timescale integrator under MPI.
- **Algorithms:** nested multiple-timescale integrator (gauge fine / Higgs coarse),
  force-gradient; frozen-length (λ→∞) RATTLE constrained scalar.
- **Representations:** general irrep engine ✅ done (`rep_general`) with per-vector fast
  path. A cached `exp(i ω·T_R)` rep_matrix path (needs a Hermitian eigensolver for log U)
  would help only for pathological reps with N^n ≫ d²; for typical irreps the per-vector
  apply is already competitive.
- **Observables/IO:** susceptibility, Binder cumulant, Polyakov, FMS correlators,
  monopoles; checkpoint + ILDG/LIME.
- **Benchmarks:** reproduce the pure-SU(2) plaquette table, the SU(2)-Higgs endpoint
  (arXiv:0911.1721 / hep-lat/9809122), and the adjoint Georgi–Glashow diagram.
