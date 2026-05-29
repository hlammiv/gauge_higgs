# Conventions — Single Source of Truth

All factor/sign choices are fixed here and used verbatim throughout the code. If a
formula in a paper disagrees, convert it to *these* conventions before comparing.
Derived in full in `theory_notes.md` §0–§3.

## Algebra
- Spacetime dimension `D ≥ 2`. Directions `μ = 0..D−1`. Lattice spacing `a = 1`.
- Gauge group `SU(N)`. Generators `Tᵃ` (a = 0..N²−2) **Hermitian**, `Tr(Tᵃ Tᵇ) = ½ δᵃᵇ`.
  Built as generalized Gell-Mann matrices / 2 (Bertlmann–Krammer ordering:
  symmetric pairs, antisymmetric pairs, then diagonal).
- Links `U_μ(x) = exp(i ωᵃ Tᵃ) ∈ SU(N)`, fundamental, N×N complex. `U_{−μ}(x) = U_μ(x−μ̂)†`.
- Momenta: gauge `P_μ(x) ∈ su(N)` stored as **real components** `Pᵃ` (P = PᵃTᵃ, Hermitian);
  scalar `π_x` real, same shape as φ. Kinetic terms `½ Σ (Pᵃ)²` and `½ π·π`.
- TA projection: `[M]_TA = ½(M − M†) − (1/2N) Tr(M − M†) I`.

## Action
- Wilson gauge: `S_g = (β/N) Σ_x Σ_{μ<ν} Re Tr(1 − U_{μν})`, with `β = 2N/g₀²`.
  Local form `S_g = −(β/N) Re Tr{U_μ(x) Σ_μ(x)} + const`; `Σ` = sum of `2(D−1)` staples.
- Scalar (irrep R): `S_H = Σ_x[ φ†φ + λ(φ†φ−1)² ] − κ Σ_{x,μ} 2 Re[ φ_x† D^(R)(U_μ(x)) φ_{x+μ̂} ]`.
- Bare mass: `a²m₀² = (1−2λ)/κ − 2D`; free critical hopping `κ_c = 1/(2D)`.

## Forces (component form is the implementation ground truth; FD-tested)
- `Fᵃ_g(x) = (β/N) Im Tr{ Tᵃ U_μ(x) Σ_μ(x) }`  (EOM `dPᵃ/dτ = −Fᵃ`, `U̇ = i(ΣPᵃTᵃ)U`).
- Scalar force: `−∂S_H/∂φ_x† = −φ_x[1+2λ(φ†φ−1)] + κ Σ_μ[ D^(R)(U_μ(x))φ_{x+μ̂} + D^(R)(U_μ(x−μ̂))†φ_{x−μ̂} ]`.
- Matter-staple link force: `Fᵃ_H(x) = 2κ Re[ φ_x† (i Tᵃ_R) D^(R)(U_μ(x)) φ_{x+μ̂} ]`.
- Total link force `F_g + F_H`; both projected with the same `algebra_project`. (HiRep principle.)

## Integrators
- Leapfrog `[I_P(ε/2) I_U(ε) I_P(ε/2)]^{N_md}`. Omelyan 2MN with `λ = 0.1931833`.
- Targets: acceptance 0.65–0.90; diagnostics `⟨exp(−ΔH)⟩ = 1`, reversibility ~1e-10.

## Observables
- Plaquette `P = ⟨(1/N) Re Tr U_pl⟩`. `L_φ = ⟨φ†φ⟩`. `L_link = ⟨Re φ† D^(R)(U) φ'⟩` (transition locator).
- Susceptibility `χ = V(⟨O²⟩−⟨O⟩²)`; Binder `U4 = 1 − ⟨O⁴⟩/(3⟨O²⟩²)` (→ 2/3 crossover).
- **Never** use gauge-variant `⟨φ⟩` as an order parameter (Elitzur).
