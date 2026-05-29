# Annotated Bibliography — Gauge + Higgs HMC (arbitrary D, arbitrary irrep)

Entries are grouped by theme. `[VERIFIED]`/`[UNVERIFIED]` tags are reproduced exactly as supplied. arXiv links are given only where an id exists; pre-arXiv works are marked "no arXiv id".

---

## A. HMC algorithm, integrators, tuning, autocorrelations

- **[UNVERIFIED]** S. Duane, A. D. Kennedy, B. J. Pendleton, D. Roweth, *Hybrid Monte Carlo*, Phys. Lett. B **195** (1987) 216. (no arXiv id) — **Essential.** The original HMC: guided MD + Metropolis = exact for any step size. Read once for the idea.
- **[VERIFIED]** R. M. Neal, *MCMC using Hamiltonian dynamics*, Handbook of MCMC (2011), Ch. 5 — https://arxiv.org/abs/1206.1901 — **Supporting (but read first).** The clearest pedagogical treatment: detailed balance via reversibility + volume preservation, leapfrog, energy error, tuning, diagnostics.
- **[VERIFIED]** M. Lüscher, *Computational Strategies in Lattice QCD*, Les Houches 2009 — https://arxiv.org/abs/1002.4232 — **Essential.** Authoritative su(N) HMC: algebra conventions, `H=½(π,π)+S`, EOM `U̇=πU`, Wilson gauge force `F₀^a=−ReTr{T^a U M}`, leapfrog, Metropolis. Supplies the requested gauge formulas directly.
- **[VERIFIED]** T. Takaishi, P. de Forcrand, *Testing and tuning symplectic integrators…*, Phys. Rev. E **73** (2006) 036706 — https://arxiv.org/abs/hep-lat/0505020 — **Essential.** Practical leapfrog vs minimum-norm (Omelyan) comparison; ~10× smaller energy error, ~50% speedup; integrate-positions-first slightly better.
- **[VERIFIED]** I. P. Omelyan, I. M. Mryglod, R. Folk, *Optimized Forest-Ruth- and Suzuki-like algorithms…* — https://arxiv.org/abs/cond-mat/0110585 — **Essential.** Source of the 2MN/4MN minimum-norm integrators and λ≈0.1931833.
- **[UNVERIFIED]** J. C. Sexton, D. H. Weingarten, *Hamiltonian evolution for the hybrid Monte Carlo algorithm*, Nucl. Phys. B **380** (1992) 665. (no arXiv id) — **Essential.** Reversible nested multiple-time-scale leapfrog; the basis for splitting gauge / Higgs / fermion forces.
- **[VERIFIED]** M. A. Clark, A. D. Kennedy, *Speeding up HMC with better integrators*, PoS LAT2007:038 — https://arxiv.org/abs/0710.3611 — **Important.** Shadow ("nearby") Hamiltonian via Poisson brackets; principled integrator tuning.
- **[VERIFIED]** A. D. Kennedy, M. A. Clark, P. J. Silva, *Force Gradient Integrators*, PoS LAT2009:021 — https://arxiv.org/abs/0910.2950 — **Important.** Adds `{S,{S,T}}` term for effective 4th-order energy conservation with one extra force-gradient evaluation.
- **[VERIFIED]** M. A. Clark, B. Joo, A. D. Kennedy, P. J. Silva, *Improving dynamical lattice QCD … Poisson brackets and a force-gradient integrator*, Phys. Rev. D **84** (2011) 071502(R) — https://arxiv.org/abs/1108.1828 — **Important.** Full shadow-Hamiltonian tuning + force-gradient integrator.
- **[VERIFIED]** H. B. Meyer et al., *Exploring the HMC trajectory-length dependence of autocorrelation times…*, CPC **176** (2007) 91 — https://arxiv.org/abs/hep-lat/0606004 — **Important.** Longer-than-conventional trajectories cut autocorrelations ~2×.
- **[VERIFIED]** A. Beskos, N. Pillai, G. Roberts, J. M. Sanz-Serna, A. Stuart, *Optimal tuning of the hybrid Monte Carlo algorithm*, Bernoulli **19** (2013) — https://arxiv.org/abs/1001.4460 — **Important.** Rigorous asymptotic optimal acceptance 0.651 and `dt ∝ d^{−1/4}`.
- **[VERIFIED]** S. Schaefer, R. Sommer, F. Virotta, *Critical slowing down and error analysis in lattice QCD simulations*, Nucl. Phys. B **845** (2011) 93 — https://arxiv.org/abs/1009.5228 — **Essential.** Topological-charge freezing (z≈5 for Q²); slow-mode-aware error analysis.
- **[VERIFIED]** M. Lüscher, S. Schaefer, *Topology, the Wilson flow and the HMC algorithm*, PoS Lattice2010:015 — https://arxiv.org/abs/1009.5877 — **Supporting.** HMC trajectory dynamics and topology barriers.
- **[UNVERIFIED]** C. Liu, A. Jaster, K. Jansen, *Instabilities in Molecular Dynamics Integrators… HMC* — (cited as hep-lat/0110047) — **Supporting.** Leapfrog step-size instability threshold / acceptance collapse.
- **[UNVERIFIED]** C. Urbach, K. Jansen, A. Shindler, U. Wenger, *A Novel Multiple-Time Scale Integrator…* — (cited as arXiv:1101.0651) — **Supporting.** Generalizes Sexton–Weingarten to flexible ratios.
- **[UNVERIFIED]** Clark, Kennedy et al., *Tuning HMC using molecular dynamics force variances*, CPC (2018) — (cited as arXiv:1801.06412) — **Supporting.** Balance per-force MD-force variances to set timescales.
- **[UNVERIFIED]** E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, 2nd ed., Springer (2006). (no arXiv id) — **Supporting.** Definitive on symplectic/splitting integrators and the modified (shadow) Hamiltonian.
- **[UNVERIFIED]** S. Schaefer, *Simulations with the HMC algorithm: implementation and data analysis* (Les Houches). (no arXiv id) — **Supporting.** Step-by-step HMC implementation tutorial (worked for φ⁴; substitute SU(N) force/exp).

## B. Gauge sector: Wilson action, su(N) exponentiation, link updates

- **[VERIFIED]** C. Morningstar, M. Peardon, *Analytic smearing of SU(3) link variables*, Phys. Rev. D **69** (2004) 054501 — https://arxiv.org/abs/hep-lat/0311018 — **Essential.** Exact SU(3) `exp(iQ)=f₀+f₁Q+f₂Q²` via Cayley–Hamilton; the TA projection. Reused for stout smearing.
- **[UNVERIFIED]** A. D. Kennedy, B. J. Pendleton, *Improved heatbath method…*, Phys. Lett. B **156** (1985) 393. (no arXiv id) — **Supporting/Important.** SU(2) accept/reject heatbath (with Fabricius–Haan); the alternative pure-gauge/angular-Higgs update.
- **[UNVERIFIED]** K. Fabricius, O. Haan, *Heat bath method…*, Phys. Lett. B **143** (1984) 459. (no arXiv id) — **Supporting.** Auxiliary-variable SU(2) heatbath; FHKP link update.
- **[UNVERIFIED]** N. Cabibbo, E. Marinari, *A new method for updating SU(N) matrices…*, Phys. Lett. B **119** (1982) 387. (no arXiv id) — **Supporting.** SU(2)-subgroup updating for SU(N) links.
- **[UNVERIFIED]** S. L. Adler, *Overrelaxation method…*, Phys. Rev. D **23** (1981) 2901. (no arXiv id) — **Supporting.** Foundational overrelaxation.
- **[UNVERIFIED]** F. Knechtli, M. Günther, M. Peardon, *Lattice QCD: Practical Essentials*, SpringerBriefs (2017). (no arXiv id) — **Important.** Concise gauge force = TA(U·staple), symplectic integrators.
- **[UNVERIFIED]** T. DeGrand, C. DeTar, *Lattice Methods for QCD*, World Scientific (2006). (no arXiv id) — **Important.** Action construction, HMC forces, link exponentiation, representations.
- **[UNVERIFIED]** C. Gattringer, C. B. Lang, *QCD on the Lattice*, LNP **788**, Springer (2010). (no arXiv id) — **Essential/Supporting.** Clean gauge-HMC derivation; covariant hopping; gauge-invariant operators.

## C. Scalar/Higgs sector: action conventions, updates, SU(2)-Higgs

- **[VERIFIED]** F. Csikor, Z. Fodor, J. Hein, K. Jansen, A. Jaster, I. Montvay, *Electroweak Phase Transition… SU(2) Higgs Model* — https://arxiv.org/abs/hep-lat/9507024 — **Essential.** Canonical variable-length SU(2)-Higgs action (β,κ,λ; Φ as 2×2 matrix); heatbath+overrelaxation+scalar-reflection.
- **[VERIFIED]** Z. Fodor, J. Hein, K. Jansen, A. Jaster, I. Montvay, *Simulating the EWPT in the SU(2) Higgs Model*, Nucl. Phys. B **439** (1995) 147 — https://arxiv.org/abs/hep-lat/9409017 — **Essential.** Hopping/2×2 form, bare-mass relation `m₀²=(1−2λ)/κ−8`, `Φ=ρ·α` with `Φ†=τ₂Φ^Tτ₂`, Bunk scalar heatbath, multicanonical.
- **[VERIFIED]** C. Bonati, G. Cossu, M. D'Elia, A. Di Giacomo, *Phase diagram of the lattice SU(2) Higgs model*, Nucl. Phys. B **828** (2010) 390 — https://arxiv.org/abs/0911.1721 — **Essential.** Fixed-length action (linear in φ → heatbath+OR, no HMC); first-order endpoint near β≈2.72; corrects earlier first-order claims to crossovers.
- **[UNVERIFIED]** C. Bonati, G. Cossu, M. D'Elia, A. Di Giacomo, *On the phase diagram of the Higgs SU(2) model* (2009) — https://arxiv.org/abs/0901.4429 — **Essential.** Explicit action; defines χ=L⁴(⟨O²⟩−⟨O⟩²), Binder V4, Z2-monopole density; Binder → 0.66667 crossover.
- **[VERIFIED]** M. Deka, S. Digal, *A New Order Parameter for the Higgs Transition in SU(2)-Higgs Theory* (2022) — https://arxiv.org/abs/2206.12884 — **Supporting.** φ as 2×2 isospin matrix; `μ₀²=(1−2λ)/κ−8`; κ,λ exactly the user's convention. Cross-check parameter map.
- **[VERIFIED]** R. Ikeda, S. Kato, K.-I. Kondo, A. Shibata, *Gauge-independent transition separating confinement and Higgs phases…* (2023) — https://arxiv.org/abs/2308.13430 — **Supporting.** Frozen unit-length action `S_H=(γ/2)ΣReTr(1−Θ†UΘ)`; gauge-invariant order parameters.
- **[VERIFIED]** M. Günther, R. Hollwieser, F. Knechtli, *Constrained Hybrid Monte Carlo algorithms for gauge-Higgs models*, CPC **251** (2020) 107081 — https://arxiv.org/abs/1908.10950 — **Supporting/Important.** Concrete (RATTLE) constrained HMC; MD treatment of the Higgs and constraint handling.
- **[UNVERIFIED]** W. Langguth, I. Montvay, P. Weisz, *Monte Carlo study of the standard SU(2) Higgs model*, Nucl. Phys. B **277** (1986) 11. (no arXiv id) — **Essential/Supporting.** Defines the standard 4D SU(2)-Higgs action; historical first-order/two-state (later revised) picture.
- **[UNVERIFIED]** B. Bunk, E.-M. Ilgenfritz, J. Kripfganz, A. Schiller, *Aspects of the SU(2)-Higgs model on the lattice*, Nucl. Phys. B **403** (1993) 453. (no arXiv id) — **Supporting.** Pioneering 4D sims; the "Bunk" scalar heatbath.
- **[UNVERIFIED]** B. Bunk, *Monte Carlo methods… electroweak phase transition*, Nucl. Phys. B (Proc. Suppl.) **42** (1995) 566. (no arXiv id) — **Important.** The radial-mode reflection update that kills ρ autocorrelations.
- **[UNVERIFIED]** W. Bock, H. G. Evertz, J. Jersak, D. P. Landau, T. Neuhaus, J. L. Xu, *Search for critical points in the SU(2) Higgs model*, Phys. Rev. D **41** (1990) 2573. (no arXiv id) — **Supporting.** Variable-length phase structure and critical-point methodology.
- **[UNVERIFIED]** I. Montvay, G. Münster, *Quantum Fields on a Lattice*, CUP (1994). (no arXiv id) — **Essential.** Authoritative gauge-Higgs conventions, hopping-parameter formulation, radial/angular decomposition, bare-mass relation.

## D. Arbitrary representations & group-theory machinery

- **[VERIFIED]** L. Del Debbio, A. Patella, C. Pica, *Higher representations on the lattice… SU(2) with adjoint fermions*, Phys. Rev. D **81** (2010) 094503 — https://arxiv.org/abs/0805.2058 — **Essential.** Building `U_R`, `T^a_R` from the fundamental; adjoint trace formula; two-index sym/antisym; Dynkin/Casimir tables.
- **[VERIFIED]** R. A. Bertlmann, P. Krammer, *Bloch vectors for qudits*, J. Phys. A **41** (2008) 235303 — https://arxiv.org/abs/0806.1174 — **Essential.** Generalized Gell-Mann basis for arbitrary SU(N): symmetric/antisymmetric/diagonal families with exact matrix elements.
- **[VERIFIED]** R. M. Fonseca, *GroupMath: A Mathematica package for group theory calculations*, CPC **267** (2021) 108085 — https://arxiv.org/abs/2011.01764 — **Essential.** Algorithms to port: RepMatrices (highest-weight), DimR, Casimir, DynkinIndex, ReduceRepProduct.
- **[VERIFIED]** R. M. Fonseca, *Calculating the RGEs of a SUSY model with Susyno*, CPC **183** (2012) 2298 — https://arxiv.org/abs/1106.5016 — **Important.** SU(N) RepMatrices with diagonal generators last; documents the generator-basis convention to match.
- **[VERIFIED]** N. Yamatsu, *Finite-Dimensional Lie Algebras and Their Representations for Unified Model Building* — https://arxiv.org/abs/1511.08771 — **Important.** Tables of dims, Dynkin indices, Casimirs, branching/projection matrices for cross-checking generated irreps.
- **[VERIFIED]** D. Bossion, P. Huo, *General Formulas of the Structure Constants in the su(N) Lie Algebra* — https://arxiv.org/abs/2108.07219 — **Supporting.** Closed-form index-only f^abc, d^abc in the GGM basis; analytic fast path / cross-check.
- **[UNVERIFIED]** H. E. Haber, *The eigenvalues of the quadratic Casimir operator…* (SCIPP notes). (no arXiv id) — **Important.** `C₂(R)=(M,M+2δ)` via inverse Cartan matrix; general-N general-irrep Casimir algorithm.
- **[UNVERIFIED]** H. Georgi, *Lie Algebras in Particle Physics*, 2nd ed. (1999). (no arXiv id) — **Important.** Highest-weight / ladder-operator construction; explicit SU(2) spin-j and SU(3) (p,q).
- **[UNVERIFIED]** M. A. A. van Leeuwen, A. M. Cohen, B. Lisser, *LiE*. (no arXiv id) — **Supporting.** Classic Lie-group computation package; conventions to mirror.

## E. Adjoint / higher-rep Higgs lattice actions

- **[VERIFIED]** C. Bonati, A. Franchi, A. Pelissetto, E. Vicari, *3D lattice SU(Nc) gauge theories with multiflavor adjoint scalars*, Phys. Rev. D **104** (2021) 094513 — https://arxiv.org/abs/2106.15152 — **Essential.** Cleanest adjoint action `Tr[Φ^t U^adj Φ]`, `U^adj,ab=2Tr(U†T^aUT^b)`; flavor order parameter `Q^{fg}`.
- **[VERIFIED]** K. Kajantie, M. Laine, A. Rajantie, K. Rummukainen, M. Tsypin, *Phase diagram of 3D SU(3)+adjoint Higgs*, JHEP **04** (1998) 023 — https://arxiv.org/abs/hep-lat/9811004 — **Essential.** Adjoint hopping `Tr[A₀ U A₀ U†]`; gauge-invariant `Tr A₀³` order parameter (N≥3); SU(2) vs SU(3) contrast.
- **[VERIFIED]** G. Catumba et al., *Lattice study of SU(2) gauge theory coupled to four adjoint Higgs fields* (2024) — https://arxiv.org/abs/2407.15422 — **Important.** Explicit `−2κ Σ Tr(Φ U Φ U†)`, `Φ=Φ^α σ^α/2`; multiple broken phases; HMC reference.
- **[VERIFIED]** A. Hart, O. Philipsen, J. D. Stack, M. Teper, *Phase diagram of the SU(2) adjoint Higgs model in 2+1 D*, Phys. Lett. B **396** (1997) 217 — https://arxiv.org/abs/hep-lat/9612021 — **Essential.** Georgi–Glashow lattice: confinement vs Higgs, photon mass, string tension, monopole density, first-order + endpoint.
- **[UNVERIFIED]** R. C. Brower, D. A. Kessler, T. Schalk, H. Levine, M. Nauenberg, *SU(2) adjoint Higgs model*, Phys. Rev. D **25** (1982) 3319. (no arXiv id) — **Essential.** Original lattice adjoint Higgs; β_H interpolates SU(2)↔U(1); two-phase structure.
- **[UNVERIFIED]** H. Georgi, S. L. Glashow, *Unified weak and EM interactions without neutral currents*, Phys. Rev. Lett. **28** (1972) 1494. (no arXiv id) — **Important.** Continuum Georgi–Glashow model.
- **[UNVERIFIED]** G. 't Hooft, *Magnetic monopoles in unified gauge theories*, Nucl. Phys. B **79** (1974) 276. (no arXiv id) — **Essential.** 't Hooft–Polyakov monopole.
- **[UNVERIFIED]** A. M. Polyakov, *Particle spectrum…* (1974) & *Quark confinement and topology…*, Nucl. Phys. B **120** (1977) 429. (no arXiv id) — **Essential.** 3D confinement by monopole condensation.
- **[UNVERIFIED]** Ph. de Forcrand, O. Jahn, *SO(3) versus SU(2) lattice gauge theory* — https://arxiv.org/abs/hep-lat/0205026 — **Supporting.** Adjoint (SO(3)) dynamics, monopoles/vortices.

## F. Phase structure, Fradkin–Shenker, observables, FMS

- **[UNVERIFIED]** E. Fradkin, S. H. Shenker, *Phase diagrams of lattice gauge theories with Higgs fields*, Phys. Rev. D **19** (1979) 3682. (no arXiv id) — **Essential.** Fundamental Higgs ⇒ Higgs and confinement analytically connected; phase boundary CAN exist for adjoint/charge-N>1.
- **[UNVERIFIED]** K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. **110** (1978) 440. (no arXiv id) — **Essential/Important.** Rigorous analyticity result underlying Fradkin–Shenker.
- **[UNVERIFIED]** J. Fröhlich, G. Morchio, F. Strocchi, *Higgs phenomenon without symmetry breaking order parameter*, Nucl. Phys. B **190** (1981) 553. (no arXiv id) — **Essential.** FMS mechanism: gauge-invariant composite operators (Tr Φ†Φ for Higgs, Tr σ^aΦ†UΦ for W) give the physical spectrum.
- **[VERIFIED]** K. Kajantie, M. Laine, K. Rummukainen, M. Shaposhnikov, *The Electroweak Phase Transition: A Non-Perturbative Analysis*, Nucl. Phys. B **466** (1996) 189 — https://arxiv.org/abs/hep-lat/9510020 — **Essential.** 3D effective SU(2)+Higgs; x,y parameters; endpoint `x_c=0.0983(15)` (m_H,c≈72 GeV).
- **[VERIFIED]** K. Kajantie, M. Laine, K. Rummukainen, M. Shaposhnikov, *Generic rules for high-T dimensional reduction…*, Nucl. Phys. B **458** (1996) 90 — https://arxiv.org/abs/hep-ph/9508379 — **Important.** Dimensional-reduction rules mapping hot 4D SM to 3D effective theory.
- **[VERIFIED]** M. Laine, A. Rajantie, *Lattice-continuum relations for 3d SU(N)+Higgs theories*, Nucl. Phys. B **513** (1998) 471 — https://arxiv.org/abs/hep-lat/9705003 — **Important.** Explicit 3D lattice action (β_G,β_H,β_R); 2-loop continuum-exact mass counterterms (fundamental and adjoint).
- **[VERIFIED]** K. Kajantie, K. Rummukainen, M. Shaposhnikov, *A Lattice MC Study of the Hot Electroweak Phase Transition*, Nucl. Phys. B **407** (1993) 356 — https://arxiv.org/abs/hep-ph/9305345 — **Supporting.** Early 3D SU(2)+adjoint+fundamental Higgs study.
- **[VERIFIED]** F. Csikor, Z. Fodor, J. Heitger et al., *Four-dimensional Simulation of the Hot EWPT with the SU(2) Gauge-Higgs Model* — https://arxiv.org/abs/hep-lat/9612023 — **Important.** Full 4D heatbath+OR production workflow.
- **[VERIFIED]** K. Kajantie, M. Laine, K. Rummukainen, M. Shaposhnikov, *Results from 3D Electroweak Phase Transition Simulations* — https://arxiv.org/abs/hep-lat/9509086 — **Important.** Heatbath+OR on the 3D effective theory; continuum/infinite-volume extrapolations.
- **[UNVERIFIED]** K. Kajantie, M. Laine, K. Rummukainen, M. Shaposhnikov, *Is there a hot electroweak phase transition at m_H ≳ m_W?*, Phys. Rev. Lett. **77** (1996) 2887 — https://arxiv.org/abs/hep-ph/9605288 — **Important.** Concise endpoint/critical-Higgs-mass statement.
- **[UNVERIFIED]** Y. Aoki, F. Csikor, Z. Fodor, A. Ukawa, *Endpoint of the first-order PT of the SU(2) gauge-Higgs model on a 4D isotropic lattice*, Phys. Rev. D **60** (1999) 013001 — https://arxiv.org/abs/hep-lat/9809122 — **Important.** Direct 4D endpoint: λ_c=0.00102(3)–0.00116(16), m_H,c=68.2(6.6) GeV.
- **[UNVERIFIED]** F. Csikor et al., *Endpoint… (Lee-Yang/Binder analysis)* — https://arxiv.org/abs/hep-lat/9901021 — **Supporting.** Companion first-order-vs-crossover methodology.
- **[UNVERIFIED]** C. Bonati, A. Pelissetto, E. Vicari, *Multicritical point of the 3D Z₂ gauge-Higgs model*, Phys. Rev. B **105** (2022) 165138 — https://arxiv.org/abs/2112.01824 — **Important.** Cheapest validation model: self-dual multicritical XY/O(2) point.
- **[UNVERIFIED]** O. Philipsen, M. Teper, H. Wittig, *On the mass spectrum of the SU(2) Higgs model in 2+1 D*, Nucl. Phys. B **469** (1996) 445 — https://arxiv.org/abs/hep-lat/9602006 — **Important.** Extracting 0⁺⁺ and 1⁻⁻ masses from smeared gauge-invariant correlators.
- **[UNVERIFIED]** T. A. DeGrand, D. Toussaint, *Topological excitations and MC simulation of Abelian gauge theory*, Phys. Rev. D **22** (1980) 2478. (no arXiv id) — **Supporting.** Standard lattice monopole definition.
- **[UNVERIFIED]** J. Greensite, K. Matsuyama, *A confinement criterion for gauge theories with matter fields* (PRD 96/98/101). (no arXiv id) — **Supporting.** Custodial/separation-of-charge criterion distinguishing confinement vs Higgs despite analytic continuity.
- **[UNVERIFIED]** L. Del Debbio, M. Faber, J. Greensite, S. Olejnik, *Center vortices and confinement vs. screening*, Phys. Rev. D **57** (1998) 2603 — https://arxiv.org/abs/hep-th/9712248 — **Important.** N-ality vs Casimir scaling; adjoint sources screened/string-broken.

## G. Combined gauge+Higgs HMC implementations (modern)

- **[VERIFIED]** B. H. Wellegehausen, A. Wipf, C. Wozar, *Phase diagram of the lattice G(2) Higgs Model*, Phys. Rev. D **83** (2011) 114502 — https://arxiv.org/abs/1102.1900 — **Essential.** Local HMC for gauge+Higgs; T≈0.75, dt≈0.25, >99% acceptance; why LHMC beats heatbath at large κ.
- **[VERIFIED]** A. D. Kennedy, K. M. Bitar, *An Exact Local Hybrid Monte Carlo Algorithm for Gauge Theories*, Nucl. Phys. B Proc. Suppl. **34** (1994) 786 — https://arxiv.org/abs/hep-lat/9311017 — **Important.** LHMC ≡ tunable overrelaxation.
- **[VERIFIED]** P. Marenzoni, L. Pugnetti, P. Rossi, *Measure of Autocorrelation Times of LHMC for Lattice QCD*, Phys. Lett. B **315** (1993) 152 — https://arxiv.org/abs/hep-lat/9306013 — **Supporting.** LHMC autocorrelation advantage.
- **[VERIFIED]** G. Catumba et al., *Lattice investigation of custodial two-Higgs-doublet model at weak quartic couplings*, JHEP **10** (2025) 214 — https://arxiv.org/abs/2507.07759 — **Important.** Modern combined gauge+scalar HMC; quaternion-parametrized doublets.
- **[VERIFIED]** G. Catumba et al., *Progress in lattice simulations for two Higgs doublet models*, PoS(LATTICE2024)145 — https://arxiv.org/abs/2412.13896 — **Supporting.** HMC setup, momenta, integrator, quaternion scalar formalism.

## H. Software architecture & performance

- **[VERIFIED]** P. A. Boyle, G. Cossu, A. Yamaguchi, A. Portelli, *Grid: A next generation data parallel C++ QCD library*, PoS LATTICE2015:023 — https://arxiv.org/abs/1512.03487 — **Essential.** Virtual-node SIMD interleaving; vReal/iMatrix/Lattice; ~350-line decltype expression-template engine; Cshift/Stencil halos.
- **[VERIFIED]** V. Drach, S. Martins, C. Pica, A. Rago, *High-Performance Simulations of Higher Representations of Wilson Fermions (HiRep v2)* — https://arxiv.org/abs/2503.06721 — **Essential.** GPU (R)HMC; SU(N) + higher reps; representation-agnostic-force architecture.
- **[UNVERIFIED]** HiRep source repository (autosun generator; MkFlags NG/REPR/GAUGE_GROUP), github.com/claudiopica/HiRep. (no arXiv id) — **Essential.** Build-time C++ generator emitting per-(N,R) headers; compile-time NG/REPR; ranlxd RNG; ILDG/openQCD/MILC converters.
- **[UNVERIFIED]** HILA framework (hilapp Clang source-to-source), github.com/CFT-HY/HILA. (no arXiv id) — **Essential.** Arbitrary compile-time NDIM; AoS/AoSoA/SoA layout switching; Parity + X iterator; adjoint/sym/antisym rep classes; gauge+complex-scalar Higgs app.
- **[VERIFIED]** L. Mazur, D. Bollweg, D. A. Clarke et al. (HotQCD), *SIMULATeQCD: A simple multi-GPU lattice code for QCD*, CPC **300** (2024) 109164 — https://arxiv.org/abs/2306.01098 — **Important.** Functor + RunFunctors kernels; centralized MemoryManagement; even-odd memory ordering; GPU-aware MPI halos.
- **[VERIFIED]** R. G. Edwards, B. Joo, *The Chroma Software System for Lattice QCD*, Nucl. Phys. B Proc. Suppl. **140** (2005) 832 — https://arxiv.org/abs/hep-lat/0409003 — **Important.** Layered QDP++ data-parallel architecture; Subsets/rb; QMP/QIO/ILDG; PETE verbosity cautionary tale.
- **[VERIFIED]** X.-Y. Jin, J. C. Osborn, *QEX: a framework for lattice field theories*, PoS LATTICE2016:271 — https://arxiv.org/abs/1612.02750 — **Supporting.** Nim metaprogramming alternative to C++ templates.
- **[UNVERIFIED]** A. Bazavov et al. (MILC), *The MILC Code*. (no arXiv id) — **Supporting.** Long-lived 4D SU(3) MIMD reference; AoS "site" struct (what to avoid for SIMD/GPU).
- **[UNVERIFIED]** M. Lüscher, *A portable high-quality RNG for lattice field theory (RANLUX)*, CPC **79** (1994) 100 — https://arxiv.org/abs/hep-lat/9309020 — **Supporting.** RANLUX (used by HiRep/openQCD); compare to counter-based.
- **[UNVERIFIED]** J. K. Salmon, M. A. Moraes, R. O. Dror, D. E. Shaw, *Parallel Random Numbers: As Easy as 1,2,3 (Random123/Philox)*, SC'11. (no arXiv id) — **Supporting.** Counter-based RNG for layout-independent reproducible per-site streams. Recommended.
