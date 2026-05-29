# Theory Notes — Implementation-Grade Reference

Self-contained notes to implement the arbitrary-D, arbitrary-N, arbitrary-irrep gauge+Higgs HMC. **Conventions are fixed once here and used throughout.**

## 0. Conventions (fix these and never deviate)

- Spacetime dimension `D` (any ≥2). Directions `μ = 0..D−1`; unit vectors `μ̂`. Lattice spacing `a=1`.
- Gauge group `G = SU(N)`. Generators `T^a` (a=1..N²−1) **Hermitian**, normalized `Tr(T^a T^b) = ½ δ^{ab}`.
- Links `U_μ(x) = exp(i ω^a T^a) ∈ SU(N)`, fundamental rep, stored as N×N complex matrices. Backward link `U_{−μ}(x) = U_μ(x−μ̂)†`.
- Wilson coupling `β = 2N/g₀²` for SU(N) (for SU(2), `β = 4/g²`).
- Scalar `φ_x` in irrep R (column vector for complex R, real vector for real R, or Hermitian-traceless matrix Φ for the adjoint).
- Momenta: gauge `P_μ(x) ∈ su(N)` (Hermitian, `P = P^a T^a`); scalar `π_x` real, same shape as φ.
- TA (traceless anti-Hermitian) projection: `[M]_TA = ½(M − M†) − (1/2N) Tr(M − M†) I`.

> **Alternative (Lüscher) anti-Hermitian convention:** `X† = −X`, `Tr(T^aT^b)=−½δ^{ab}`, `U=exp(X)`, `U̇=πU`. Physics identical; differs by factors of i. Pick the Hermitian one above for clarity and be rigorous.

## 1. The action

### 1.1 Wilson gauge action
```
S_g = (β/N) Σ_x Σ_{μ<ν} Re Tr( 1 − U_{μν}(x) ),
U_{μν}(x) = U_μ(x) U_ν(x+μ̂) U_μ(x+ν̂)† U_ν(x)†.
```
Single-link local form (isolating `U_μ(x)`):
```
S_g = −(β/N) Re Tr{ U_μ(x) Σ_μ(x) } + const,
Σ_μ(x) = Σ_{ν≠μ} [ U_ν(x+μ̂) U_μ(x+ν̂)† U_ν(x)†          (forward staple)
                  + U_ν(x+μ̂−ν̂)† U_μ(x−ν̂)† U_ν(x−ν̂) ]  (backward staple).
```
There are `2(D−1)` staples per link.

### 1.2 Scalar action (general irrep R)
Gauge-covariant hopping + on-site potential:
```
S_H = Σ_x [ φ_x†φ_x + λ(φ_x†φ_x − 1)² ]
     − κ Σ_x Σ_{μ=0}^{D−1} 2 Re[ φ_x† D^(R)(U_μ(x)) φ_{x+μ̂} ].
```
- `D^(R)(U) = exp(i ω^a T^a_R)` — the SAME group element promoted to rep R via the rep-R generators `T^a_R`. **No separate links per rep.**
- Bare-mass relation (to convert κ,λ ↔ m₀²): `a² m₀² = (1−2λ)/κ − 2D`. Free critical hopping `κ_c = 1/(2D)`.
- Field normalization: `φ_cont = √(2κ) φ_lat`.
- **Fundamental rep:** `D^(R)(U)=U`. **Real rep** (e.g. adjoint, G2-7): φ real, `D^(R)(U)` real orthogonal.
- **λ→∞ (frozen length):** `|φ|=1`; the scalar becomes a sphere/group element. Use constrained MD (RATTLE) or local heatbath.

### 1.3 SU(2)-fundamental matrix form (the canonical benchmark)
Package the doublet as `Φ = φ₀ I + i φ_a τ_a` (τ = Pauli), obeying `Φ† = τ₂ Φ^T τ₂`, so `½Tr(Φ†Φ) = φ†φ = ρ²`, factor `Φ = ρ·α`, α∈SU(2):
```
S = β Σ_pl [1 − ½ Re Tr U_pl]
  + Σ_x { ½Tr(Φ_x†Φ_x) + λ(½Tr(Φ_x†Φ_x) − 1)²
          − κ Σ_μ Tr(Φ_{x+μ̂}† U_μ(x) Φ_x) }.
```
(The single κ here = 2κ of the doublet form because Tr(Φ†Φ)=2φ†φ.)

### 1.4 Adjoint scalar (two equivalent forms)
Hermitian-traceless matrix `Φ = φ^a T^a`, conjugation transform `Φ → Ω Φ Ω†`:
```
S_H^adj ⊃ −2κ Σ_{x,μ} Tr[ Φ_x U_μ(x) Φ_{x+μ̂} U_μ(x)† ].
```
Equivalent real-vector form with the adjoint link matrix:
```
[U^adj_μ(x)]^{ab} = 2 Tr[ T^a U_μ(x) T^b U_μ(x)† ]   (real orthogonal, det +1)
S_H^adj ⊃ −κ Σ_{x,μ} φ^a_x [U^adj_μ(x)]^{ab} φ^b_{x+μ̂}  (×2 for the Hermitian-matrix normalization).
```

## 2. HMC Hamiltonian and equations of motion

Augmented (fictitious) Hamiltonian:
```
H[P,π,U,φ] = ½ Σ_{x,μ,a} (P^a_μ(x))² + ½ Σ_x π_x·π_x + S_g[U] + S_H[φ,U].
```
Target joint density `exp(−H)`; marginal over momenta reproduces `exp(−S)`.

Momentum heatbath each trajectory: `P^a ~ N(0,1)` i.i.d. (assemble `P=P^aT^a`); `π ~ N(0,1)` per scalar component.

Hamilton's equations (MD):
```
dU_μ(x)/dτ = i P_μ(x) U_μ(x),      dP_μ(x)/dτ = −(F_g + F_H)_μ(x),
dφ_x/dτ    = π_x,                  dπ_x/dτ    = −∂S_H/∂φ_x.
```
Link drift integrated as `U ← exp(i P dt) U`; scalar drift as `φ ← φ + dt·π`.

Metropolis: `P_acc = min(1, exp(−ΔH))`, `ΔH = H_final − H_initial`. Exactness needs reversible + area-preserving MD. Diagnostic: `⟨exp(−ΔH)⟩ = 1` ⇒ `⟨ΔH⟩ ≥ 0`.

## 3. Forces (all live in su(N) or the scalar tangent space)

### 3.1 Gauge force (per link, su(N)-valued)
```
F_g,μ(x) = (β/N) [ U_μ(x) Σ_μ(x) ]_TA.
```
Component form (Lüscher): `F^a = −(β/N) Re Tr{ T^a U_μ(x) Σ_μ(x) }`. Sign follows `dP/dτ = −F`.

### 3.2 Scalar force (on φ; the MD force on the scalar field)
From `S_H` (general R):
```
−∂S_H/∂φ_x† = −φ_x[1 + 2λ(φ_x†φ_x − 1)]
              + κ Σ_μ [ D^(R)(U_μ(x)) φ_{x+μ̂} + D^(R)(U_μ(x−μ̂))† φ_{x−μ̂} ].
```
The bracketed sum is the gauge-covariant "hopping staple" over the 2D neighbors. For frozen-length scalars, project the force onto the constraint tangent space (or use RATTLE).

### 3.3 Scalar back-reaction on the gauge link (the "matter staple")
Differentiate the hopping term w.r.t. the fundamental link using `∂/∂ω^a D^(R)(U) = i T^a_R D^(R)(U)`:
```
F_H,μ(x)^a = −∂S_H/∂ω^a = 2κ Re[ φ_x† (i T^a_R) D^(R)(U_μ(x)) φ_{x+μ̂} ].
```
- **Fundamental:** matrix form `F_H,μ(x) = 2κ [ U_μ(x) φ_{x+μ̂} φ_x† ]_TA` (a rank-1 matter staple added to the gauge force).
  - SU(2) matrix convention: `F_H,μ(x) = κ [ U_μ(x) Φ_x Φ_{x+μ̂}† ]_TA`.
- **Adjoint:** reduces to the TA projection of the scalar staple `Σ_μ^Φ = Φ_{x+μ̂} ... ` (sum over flavors if multiflavor).
- **General R:** assemble `F^a` via the `iT^a_R` kernel above, then store as an su(N) algebra vector (the same `_algebra_project` routine used by the gauge force).
- **Robust fallback:** reverse-mode autodiff of `S_H` reproduces the force for arbitrary R with no hand derivation.

The total link force is `F_g + F_H`, both ending in the same TA/algebra projection — this is the HiRep representation-agnostic principle.

## 4. Symplectic integrators

Elementary maps: `I_P(ε): P ← P − εF` (kick), `I_U(ε): U ← exp(iεP)U` and `φ ← φ + εint` (drift).

### 4.1 Leapfrog (Störmer–Verlet, 2LF)
```
J(ε) = [ I_P(ε/2) I_U(ε) I_P(ε/2) ]^{N_md},   ε = τ/N_md.
```
Reversible, symplectic; `ΔH = O(ε²)` (no secular drift; a shadow Hamiltonian is exactly conserved).

### 4.2 2nd-order minimum-norm (Omelyan/2MN), position-update form
```
U_2MN(ε) = e^{λε S} e^{(ε/2)T} e^{(1−2λ)ε S} e^{(ε/2)T} e^{λε S},
λ ≈ 0.1931833   (minimum-norm; lattice-optimal value mildly parameter-dependent).
```
S = force kick, T = drift. ~10× smaller RMS energy error than leapfrog ⇒ ~3× larger ε ⇒ ~50% net speedup.

### 4.3 Nested multiple-timescale (Sexton–Weingarten)
Put cheap/stiff force on fine substeps, expensive/soft on coarse:
```
U(ε_coarse) = [e^{(ε_c/2)S_slow}] [U_fine(ε_c/m)]^m [e^{(ε_c/2)S_slow}],   ε_fine = ε_c/m, m integer.
```
For gauge+Higgs: gauge force on the fine scale (it is stiffer), Higgs force on the coarse scale; nest fermions on a yet-coarser scale later. Tune `m` by balancing per-force MD-force variance.

### 4.4 Force-gradient (optional, effective O(ε⁴))
Add `−c·ε³{S,{S,T}}` to a kick (directional Hessian of S along F). Cancels leading shadow error with one extra force-gradient evaluation. Structure the integrator as a composition of reusable kick/drift/FG-kick operators so it drops in.

### 4.5 Tuning
Target acceptance ~0.65–0.90 (asymptotic optimum 0.651, `ε ∝ d^{−1/4}`). Scan trajectory length τ (start τ=1) to minimize integrated autocorrelation of observables; longer τ often helps ~2×.

## 5. su(N) exponentiation (link drift)

- **SU(2):** `exp(i θ n̂·σ) = cos|θ| I + i sin|θ| (n̂·σ)`; closed form, stays exactly in SU(2). (Quaternion fast path recommended.)
- **SU(3):** Morningstar–Peardon Cayley–Hamilton, `exp(iQ)=f₀I+f₁Q+f₂Q²`, with `c₀=det Q=(1/3)Tr Q³`, `c₁=½Tr Q²`, eigenvalues `2u, −u±w`, `u=√(c₁/3)cos(θ/3)`, `w=√c₁ sin(θ/3)`, `θ=arccos(c₀/c₀^max)`, `c₀^max=2(c₁/3)^{3/2}`. Use the small-w stable expansion of `ξ₀(w)=sin w/w`. No diagonalization.
- **General SU(N):** scaled Taylor series + scaling-and-squaring, or Hermitian eigendecomposition of `H=ω^aT^a` then `U=exp(iH)`.
- **Reunitarize** every trajectory (or every few): modified Gram–Schmidt on rows, then rescale to det=1; log `max‖U†U−I‖`.

## 6. Representation-matrix construction (the arbitrary-irrep engine)

### 6.1 Generators (generalized Gell-Mann, any SU(N))
Three families (Bertlmann–Krammer): symmetric `Λ^{jk}_s=|j⟩⟨k|+|k⟩⟨j|`, antisymmetric `Λ^{jk}_a=−i|j⟩⟨k|+i|k⟩⟨j|` (1≤j<k≤N), diagonal `Λ^l=√(2/(l(l+1)))(Σ_{j≤l}|j⟩⟨j| − l|l+1⟩⟨l+1|)` (1≤l≤N−1). Normalized `Tr(ΛΛ)=2δ`; set `T^a_F=Λ^a/2`. (N=2 → Pauli/2; N=3 → Gell-Mann.) Structure constants by trace: `f^{abc}=−2i Tr([T^a,T^b]T^c)`, `d^{abc}=2 Tr({T^a,T^b}T^c)`.

### 6.2 Building T^a_R and D^(R)(U)
- **Adjoint (fast path, no exp):** `D^adj(U)_{ab}=2 Tr(T^a U T^b U†)` (real, orthogonal). Generators `(T^a_adj)_{bc}=−i f^{abc}`.
- **Two-index sym/antisym:** `(T^a_{2S/2A})_{ij,kl}=(1/√2)(T^a_{F,ik}δ_{jl} ± T^a_{F,il}δ_{jk})` (+ symmetric dim N(N+1)/2, − antisym dim N(N−1)/2). The link can be formed directly as a symmetrized bilinear of U.
- **General irrep:** highest-weight/weight-space construction (Georgi; GroupMath RepMatrices): diagonalize Cartan generators to weights, build ladder operators with Chevalley normalization. Tag each irrep by Dynkin labels.
- To get `ω^a` from a fundamental U: `H=−i log U` (eigendecomp of U), `ω^a=2 Tr(T^a_F H)`; then `D^(R)(U)=exp(i ω^a T^a_R)`.

### 6.3 Invariants (runtime self-checks)
`Tr(T^a_R T^b_R)=T(R)δ^{ab}`, `T^a_R T^a_R = C₂(R)I`, `d(R)C₂(R)=T(R)(N²−1)`.
- fund: d=N, T=½, C₂=(N²−1)/2N.
- adjoint: d=N²−1, T=N, C₂=N.
- 2-sym: d=N(N+1)/2, T=(N+2)/2, C₂=(N−1)(N+2)/N.
- 2-antisym: d=N(N−1)/2, T=(N−2)/2, C₂=(N+1)(N−2)/N.
- General: `C₂(R)=(M,M+2δ)` via inverse Cartan matrix (Haber). Cross-check against Yamatsu tables.

### 6.4 N-ality (controls phase structure)
`D^(R)(z·1)=z^{k(R)}·1`, z=e^{2πi/N}, `k(R)`=N-ality. Adjoint k=0 (center-blind ⇒ genuine transition possible); fundamental k=1 (Fradkin–Shenker continuity). Make N-ality a first-class attribute of each Representation object.

## 7. Observables (validation + physics)

- Average plaquette `P = ⟨(1/N) Re Tr U_pl⟩` (SU(2): `½ Re Tr`). Strong coupling `P~β/2N`; weak `P→1`.
- Higgs length `L_φ = ⟨φ†φ⟩` (= ⟨ρ²⟩); constant in frozen limit.
- Gauge-invariant link/hopping energy `L_link = ⟨Re φ_x† D^(R)(U) φ_{x+μ̂}⟩` (SU(2): `⟨½Tr Φ†UΦ⟩`). **Most sensitive Higgs-transition locator.**
- Polyakov loop `⟨(1/N)Tr Π_t U_0⟩` (screened by matter, not a strict order parameter, but standard).
- Susceptibility `χ(O)=V(⟨O²⟩−⟨O⟩²)`; Binder cumulant `U4 = 1 − ⟨O⁴⟩/(3⟨O²⟩²)`. First order: χ peak ∝ V, U4 dips below 2/3. Crossover: χ peak saturates, U4 → 0.6667.
- Spectrum (FMS): 0⁺⁺ from `Tr Φ†Φ` / `Re Tr Φ†UΦ`; 1⁻⁻ from `O_μ^a = Tr(σ^a Φ†U_μΦ)`. Mass from exponential decay of zero-momentum connected correlators; APE/HYP smear and use a small variational basis.
- Adjoint extras: `Tr Φ²` (always), `Tr Φ³` (order parameter only for N≥3), `Tr[Φ U Φ U†]`; multiflavor `Q^{fg}=Σ_a φ^{af}φ^{ag} − (δ^{fg}/N_f)Tr`.
- Monopoles: DeGrand–Toussaint construction; Z₂-monopole density `M=1−(1/N_c)Σ_c σ_c` for bulk diagnostics.
- **Elitzur's theorem:** never use a gauge-variant `⟨φ⟩` as an order parameter.

## 8. Key physics facts to reproduce

- **Fradkin–Shenker:** fundamental matter ⇒ Higgs and confinement analytically connected (first-order line ending in a critical endpoint, else crossover). Adjoint / trivial-N-ality matter ⇒ a genuine transition can exist.
- **Benchmarks:** frozen-length 4D SU(2) endpoint near β≈2.72–2.73 (arXiv:0911.1721); 3D EWPT endpoint `x_c=0.0983(15)`, m_H,c≈72 GeV (arXiv:hep-lat/9510020); 4D isotropic endpoint λ_c≈0.001, m_H,c≈68 GeV (hep-lat/9809122); SU(2)/SU(3)-adjoint Georgi–Glashow genuine transitions (hep-lat/9612021, hep-lat/9811004).
