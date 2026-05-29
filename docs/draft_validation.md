# Master draft validation

Consolidated group-theory / linear-algebra validation of the discrete-symmetry draft
(`literature/discrete_symmetry/draft/draft-1.tex`, adjudicated where needed against
`literature/discrete_symmetry/notes_WL/notes1.tex`) for the program **break SU(N) to a
finite nonabelian crystal-like subgroup H via a single large Higgs irrep**.

Seven target breakings:

| G | H | Higgs irrep sigma | Young | dim_R | \|H\| |
|---|---|---|---|---|---|
| SU(2) | BT (2T) | spin-3 | {6} | 7 | 24 |
| SU(2) | BO (2O) | spin-4 | {8} | 9 | 48 |
| SU(2) | BI (2I) | spin-6 | {12} | 13 | 120 |
| SU(3) | Sigma(108) | (2,2) | {4,2} | 27 | 108 |
| SU(3) | Sigma(216) | (4,1) | {5,1} | 35 | 216 |
| SU(3) | Sigma(648) | (3,3) | {6,3} | 64 | 648 |
| SU(3) | Sigma(1080) | (6,0) | {6} | 28 | 1080 |

All results below are **pure group theory, linear algebra, and mean-field root-finding on
small matrices** (rep dim <= 64 for full matrix builds; characters for multiplicity
enumeration). **No HMC, no lattice simulation.** Each section corresponds to one validation
task with its own test (all passing) and its own standalone results doc.

| task | test | status | doc |
|---|---|---|---|
| Singlet multiplicities (characters) | `test_singlets` | 437 passed, 0 failed | `draft_validation_singlets.md` |
| Gauge-boson mass M_ab | `test_draft_gaugemass` | 60 passed, 0 failed | `draft_validation_gaugemass.md` |
| Scalar mass spectra | `test_draft_spectra` | 89 passed, 0 failed | `draft_validation_spectra.md` |
| Mean-field strong-weak transition | `test_meanfield` | 48 passed, 0 failed | `draft_validation_meanfield.md` |
| SU(3) bounded-from-below (BFB) | `test_bfb` | 11 passed, 0 failed | `draft_validation_bfb.md` |

---

## 1. Executive summary

### What our code CONFIRMS

1. **Singlet alignment (Observation 1).** Character projection
   `m_{sigma->tr} = (1/|H|) sum_h chi^sigma(h)` confirms that each chosen Higgs irrep is the
   **first** non-trivial irrep of its subgroup carrying an H-singlet, and that this lowest
   singlet is **unique (m = 1)** for all seven targets. The mult-2 case
   (2T, j=6, d=13, m=2) and the higher-irrep flags ((7,1)/(1,7) -> Sigma108 m=2, Sigma216 m=1;
   (4,4) -> Sigma108 m=5, Sigma216 m=3, Sigma648 m=1, Sigma1080 m=1) all reproduce exactly.
   SU(2) is center-blind: all half-integer j give m = 0. Matches `su2subduction` (j<=10) and
   `su3_table` (p+q<=6) entry-by-entry.

2. **Complete continuous-symmetry breaking + degenerate gauge bosons.** The gauge mass matrix
   `M_ab = Re(phi_B^dag T_b T_a phi_B)` at the H-singlet VEV is **proportional to the identity**
   (off-diagonal ~ 1e-14, eigenvalue spread ~ 1e-13) for all seven subgroups: every gauge boson
   is degenerate, M is rank-full over all dim_G generators, and **no continuous subgroup
   survives**. `Tr M = C2(sigma)` is reproduced as an exact identity, and the common eigenvalue
   `C2/dim_G` gives `m_A^2 = (2/dim_G) C2`. The SU(3) draft brackets
   `p^2+q^2+pq+3p+3q = 24/36/45/54` reproduce exactly. No discrepancies.

3. **Scalar spectra (Observations 2 & 3).** In every G(x)Gbar channel the eigenvalues cluster
   into H-irrep multiplets; the radial (VEV) direction is its own eigenvalue `lambda_0^(c)`; and
   the broken-generator Goldstones span the **full adjoint** (rank 3 for SU(2), rank 8 for SU(3))
   with per-channel Hessian eigenvalue **exactly `(1/3) lambda_0^(c)`** in every channel of all
   seven cases (the rigorous form of the draft's "1/3"). The full tuned complex-basis Hessian has
   exactly `dim(G)+1` zero modes (3+1 for SU(2), 8+1 for SU(3): gauge Goldstones + global U(1)).
   The **SU(2) BT/BO/BI tables reproduce the draft exactly** once the per-channel CG factor
   `1/sqrt(2J+1)` and the `lambda_0/3` subtraction are applied. The BO self-flag "check T_{1,2}"
   is **resolved**: T_1 = Goldstone, T_2 = massive.

4. **Mean-field strong-weak transition.** All eight SU(2)-sector criticals
   (beta_c, alpha_B, v_B for D=3,4; continuous SU(2) and 2T/2O/2I) reproduce the draft to
   <= 1.1e-4 (target 1e-3). The 2T closed form `(3 + 8 cosh(a/2) + cosh(a))/12` matches the group
   sum to 1e-12 and the SU(2) Bessel form `2 I_1(a)/a` matches the torus quadrature to 1e-6. The
   **discrete-frozenness signature** is confirmed: continuous G has `1 - v_B` decaying as a power
   law, discrete H exponentially (2T ~ e^{-alpha/2}; Sigma108 ~ 2.5e-8 at alpha=30).

5. **SU(3) bounded-from-below (open draft TODO `\yy`).** Resolved: BFB <=> **copositivity** of
   the channel couplings {f_rho'}, equivalently `min_{|Phi|=1} V_4 >= 0`. Confirmed:
   `all f_c = 1 -> V_4 == 1` exactly (since `sum_c P_c = 1`); the Sigma(108) locking couplings are
   BFB (`min V_4 = +0.0521`); `all f_c >= 0` is sufficient but **not** necessary (the highest-C2
   channel is load-bearing; the singlet channel is forgiving).

### Where DISCREPANCIES were found (and adjudicated)

- **`draft-1.tex` SU(3) scalar-spectrum tables (S108mass / S216mass / S648mass / S1080mass) are
  WRONG.** They disagree with our computation **and** with the project's own
  `notes_WL/notes1.tex`, failing even **convention-independent within-channel ratio** checks
  (e.g. draft S108 (3,3) `2^(0)/radial = 0.536` vs the physical `0.69048`; draft S1080 (4,4)
  `8_massive/radial = -0.016` vs the physical `0.32275`). **Our computation adjudicates in favor
  of the NOTES**, which it reproduces to 5-6 digits in every channel. The `draft-1.tex` SU(3)
  tables appear to be a stale/garbled version superseded by the notes. **Use the notes SU(3)
  tables.** (This is a draft-vs-notes table disagreement that our code resolved.)

- **`draft-1.tex` BT J=6 entries are slightly off.** The draft prints `0.60553` (A0) and
  `0.29424` (T_massive); the exact values are `24/(11 sqrt13) = 0.60513` and
  `35/(33 sqrt13) = 0.29416`, which we reproduce. The draft's last two digits are a
  rounding/transcription artifact, not a structural error.

### Non-errors worth flagging (convention points, not discrepancies)

- **Gauge-mass factor of 2.** The raw matrix diagonal is `M_aa = C2/dim_G`; the physical
  `m_A^2 = (2/dim_G) kappa g^2 C2 = 2 M_aa` carries the standard kinetic-normalization factor of
  2. Both forms are reported; the physical `m_A^2` is used for the type-I/II comparison.

- **Per-channel Hessian normalization (SU(2)).** The draft's `H^(J)` differs from our
  orthonormal-projector `V_c` Hessian by the per-channel CG singlet coefficient
  `1/sqrt(2J+1)`. After that rescale (plus the `lambda_0/3` subtraction and a single overall
  constant fixing J=0 A0 = 1) the SU(2) tables match exactly.

- **Goldstone-count basis subtlety.** The per-channel TABLE reproduction uses the draft's
  real-dof (W-basis) treatment; the full Goldstone COUNT (dim_G + 1 zero modes) requires the
  rep-faithful complex (2d-dof) treatment. Both are consistent.

- **Adjoint-Casimir channel merging.** `CasimirChannels` groups G(x)Gbar by adjoint-Casimir
  eigenvalue, merging distinct (p,q) with equal C2 (e.g. (4,1) & (1,4) at C2=12; (4,4) & (6,0) at
  C2=24). One C2-channel can therefore carry several draft (p,q) columns. Not an error.

- **BFB criterion is copositivity, not matrix positivity.** For SU(2) integer-spin (real) reps the
  quartic invariants collapse to quadratics in field bilinears (BFB trivial); for SU(3) `V_4` is
  genuinely quartic in Phi, so BFB is copositivity of {f_rho'} on the rank-1 cone.

---

## 2. Per-task key numbers

### 2.1 H-singlet multiplicity enumeration (characters only)

`m_{sigma->tr} = (1/|H|) sum_h chi^sigma(h)`, computed from traces only (no rep build). SU(2)
via `chi^j = sin((2j+1)theta)/sin(theta)`; SU(3) via the Jacobi-Trudi determinant (chosen over
the bialternant for numerical stability near degenerate eigenvalues).

Lowest non-trivial singlet of each group, all m = 1 (and each is the FIRST such irrep):

| G | H | irrep | dim | m |
|---|---|---|---|---|
| SU(2) | 2T | j=3 ({6}) | 7 | 1 |
| SU(2) | 2O | j=4 ({8}) | 9 | 1 |
| SU(2) | 2I | j=6 ({12}) | 13 | 1 |
| SU(3) | Sigma108 | (2,2) | 27 | 1 |
| SU(3) | Sigma216 | (4,1) | 35 | 1 |
| SU(3) | Sigma648 | (3,3) | 64 | 1 |
| SU(3) | Sigma1080 | (6,0) | 28 | 1 |

Multiplicity-2: 2T j=6 (d13) m=2. Higher flags: (7,1)/(1,7) d80 -> Sigma108 m=2, Sigma216 m=1;
(4,4) d125 -> Sigma108 m=5, Sigma216 m=3, Sigma648 m=1, Sigma1080 m=1. SU(2) half-integer j all
m=0 (center-blind). Matches `su2subduction` (j<=10) and `su3_table` (p+q<=6) exactly.

### 2.2 Gauge-boson mass M_ab

`M_ab = Re(phi_B^dag T_b T_a phi_B)`, H-singlet VEV unit-normalized, kappa g^2 = 1. M = (C2/dim_G) I
for all seven; `m_A^2 = (2/dim_G) C2`. SU(2): C2(j)=j(j+1), dim_G=3. SU(3):
C2(p,q)=(1/3)(p^2+q^2+pq+3p+3q), dim_G=8.

| subgroup | C2(sigma) | common eig (C2/dim_G) | m_A^2 = 2*eig | m_A |
|---|---|---|---|---|
| BT (2T, spin-3) | 12 | 4.00000 | 8.00000 | 2.82843 |
| BO (2O, spin-4) | 20 | 6.66667 | 13.33333 | 3.65148 |
| BI (2I, spin-6) | 42 | 14.00000 | 28.00000 | 5.29150 |
| Sigma108 (2,2) | 8 | 1.00000 | 2.00000 | 1.41421 |
| Sigma216 (4,1) | 12 | 1.50000 | 3.00000 | 1.73205 |
| Sigma648 (3,3) | 15 | 1.87500 | 3.75000 | 1.93649 |
| Sigma1080 (6,0) | 18 | 2.25000 | 4.50000 | 2.12132 |

SU(3) brackets `p^2+q^2+pq+3p+3q = 24/36/45/54` -> `m_A^2 = (1/12)*bracket = 2/3/3.75/4.5`,
reproduced exactly. No discrepancies.

### 2.3 Scalar mass spectra

SU(2) tables reproduce `draft-1.tex` exactly (after the `1/sqrt(2J+1)` rescale + `lambda_0/3`
subtraction + single normalization constant). Representative entries:

- **BT** J2: T_massive = sqrt5/3 = 0.74536; J4: A0 = 14/11 = 1.27273, T = -10/11 = -0.90909;
  J6: A0 = 24/(11 sqrt13) = 0.60513, T = 35/(33 sqrt13) = 0.29416 (draft's 0.60553/0.29424 are
  rounding artifacts).
- **BO** J2: T2_massive = sqrt5/11 = 0.20328, E = 4 sqrt5/11 = 0.81312; J4: A1 = 98/143,
  T2 = -10/143, E = -40/143. (T_1 = Goldstone, T_2 = massive: self-flag resolved.)
- **BI** full J0..J12; e.g. J6: A0 = 121 sqrt13/323 = 1.35069, J = -7 sqrt13/34 = -0.74232,
  H = -63 sqrt13/323 = -0.70325.

SU(3) (adjudicated to the NOTES; `draft-1.tex` SU(3) tables are WRONG), per-channel
radial-normalized (radial = 1, Goldstone = 1/3):

- **Sigma108** (3,3): 2^0 = 0.69048, massive {0.31829, 0.05076}; (4,4): 0.43333, {0.58045, 0.41538};
  (1,4): -1/6, {0.19514, -0.02847}.
- **Sigma1080** (4,4): 8m = 0.32275, 5+5' {0.25645, 0.02503}, 9 {-1/9, 0.03704}; (5,5): 8m = 0.65741;
  (6,6): 8m = 0.08838, 9 {0.33586, 0.28283}.

Lightest non-Goldstone scalar mass (representative stable locking from `find_stable_couplings`;
**coupling-dependent, not unique**):

| subgroup | zero modes | m_H^2 (lightest) |
|---|---|---|
| BT | 4 (3+U(1)) | 0.467 |
| BO | 4 (3+U(1)) | 0.125 |
| BI | 4 (3+U(1)) | 0.215 |
| Sigma108 | 9 (8+U(1)) | 0.182 |
| Sigma1080 | 9 (8+U(1)) | 0.065 |

(Sigma216 and Sigma648 were validated structure-only; the notes give no numeric table and m_H is
not a unique number, so no m_H^2 is quoted for them.)

### 2.4 Mean-field strong-weak coupling transition

Saddle point in temporal gauge: `N_P/N_l = (D-2)/2`, `E_P(v) = 1 - v^4 + (1-v^2) 2/(D-2)`,
`v_B = w'(alpha)`, `alpha = -beta (N_P/N_l) E_P'(v_B)`, `F_P = E_P - (1/beta)(N_l/N_P) s`;
beta_c is the first-order point where the nontrivial-branch F_P meets F_P^triv = E_P(0).

SU(2)-sector criticals reproduce the draft to <= 1.1e-4 (target 1e-3). D=4 examples:

| group | beta_c (ref) | alpha_B (ref) | v_B (ref) |
|---|---|---|---|
| SU(2) cont | 1.68173 (1.6817) | 4.71901 (4.7190) | 0.70430 (0.7043) |
| 2T (BT) | 1.53740 (1.5374) | 7.73094 (7.7309) | 0.92604 (0.9260) |
| 2O (BO) | 1.67377 (1.6738) | 5.05231 (5.0523) | 0.73029 (0.7303) |
| 2I (BI) | 1.68171 (1.6817) | 4.72008 (4.7201) | 0.70438 (0.7044) |

SU(3) sector (NEW, same validated machinery), D=4: SU(3) cont beta_c=5.07223; Sigma108 2.34019,
Sigma216 2.68704, Sigma648 3.23397, Sigma1080 3.48993 -- discrete groups transition at smaller
beta_c, with v_B ~ 0.99 (links lock onto group elements) vs continuous SU(3) v_B ~ 0.77. beta_c
grows monotonically with |H| (108<216<648<1080).

### 2.5 SU(3) bounded-from-below

`V_4(Phi) = sum_rho' f_rho' ||P_rho' M||^2`, M = Phi Phi^dag, homogeneous degree-4
=> BFB <=> `min_{|Phi|=1} V_4 >= 0` (copositivity). Sigma108 (2,2) channels
C2 = 0,3,6,8,12,15,18,20,24.

| coupling set | min_{|Phi|=1} V_4 | BFB |
|---|---|---|
| all f_c = 1 | +1.000000 (= \|\|M\|\|_F^2 exactly) | YES |
| f_last = -1, rest 0.05 | -0.490339 | NO |
| Sigma108 locking f | +0.052109 | YES |
| find_stable_couplings f | +0.052105 | YES |

The highest-C2 channel is load-bearing for boundedness; the singlet channel can go strongly
negative while V_4 stays bounded (provided f_last >= 0).

---

## 3. Type-I vs Type-II classification

Comparison of the degenerate gauge-boson mass `m_A` (section 2.2) to the lightest non-Goldstone
scalar mass `m_H = sqrt(m_H^2)` (section 2.3). Following the superconductor analogy:
**type I if m_A > m_H, type II if m_A < m_H**.

| subgroup | m_A | m_H = sqrt(m_H^2) | m_A / m_H | class |
|---|---|---|---|---|
| BT (2T, spin-3) | 2.82843 | 0.68337 | 4.139 | **I** |
| BO (2O, spin-4) | 3.65148 | 0.35355 | 10.328 | **I** |
| BI (2I, spin-6) | 5.29150 | 0.46368 | 11.412 | **I** |
| Sigma108 (2,2) | 1.41421 | 0.42661 | 3.315 | **I** |
| Sigma1080 (6,0) | 2.12132 | 0.25495 | 8.321 | **I** |
| Sigma216 (4,1) | 1.73205 | (not computed) | -- | -- |
| Sigma648 (3,3) | 1.93649 | (not computed) | -- | -- |

All five fully computed subgroups come out **type I** (m_A > m_H, ratio > 1) at the representative
locking couplings.

**IMPORTANT CAVEAT (this is a finding, not a firm prediction).** The ratio `m_A/m_H` above is
**convention- and coupling-dependent**, so the type-I/II verdict here is provisional:

- `m_A` is computed with a **unit-normalized** H-singlet VEV (`|phi_B| = 1`) and `kappa g^2 = 1`.
  For a physical VEV of length v, `m_A` scales by v (and the gauge coupling g enters linearly).
- `m_H^2` comes from a **representative stable-locking quartic** (`find_stable_couplings`); it is
  **not unique** -- it depends on the chosen channel couplings {f_c}, whose overall scale is not
  tied to the same VEV normalization used for `m_A`.

Because the two masses are extracted from independently-normalized inputs, the absolute ratio (and
hence the type-I/II boundary) is **not fixed at the group-theory level**. The robust, convention-
independent statement is the *qualitative ordering within each family* (larger |H| / larger C2
gives larger m_A; the scalar gap m_H is small and tunable). A firm type-I/II prediction requires
fixing v, g, and the physical {f_c} together -- a **large-numerics (lattice) task** (see section 4).
Sigma216 and Sigma648 lack a quoted m_H (notes give no numeric scalar table; m_H is coupling-
dependent), so no ratio is given for them.

---

## 4. Remaining for the lattice (large-numerics) phase

The validation above is purely tree-level / group-theoretic. The following require Monte Carlo /
HMC simulation of the full Fradkin-Shenker gauge+Higgs action (out of scope here):

1. **Fix the physical scale for type-I/II.** Determine v, g, and the physical channel couplings
   {f_c} self-consistently, then re-evaluate m_A/m_H to convert the provisional type-I verdicts
   into firm classifications. (The group-theory phase only fixes ratios up to the two independent
   normalizations.)

2. **Locate the actual phase transition.** The mean-field beta_c values (section 2.4) are tree-
   level saddle points. Confirm the order and location of the SU(N) -> H transition by full lattice
   measurement (plaquette, link, susceptibilities), and check whether the first-order mean-field
   prediction survives fluctuations.

3. **Confinement / screening observables.** Wilson and Polyakov loops, string tension,
   center-blind vs center-sensitive reps. Per the project memory, 4D discrete-H is deconfined, so
   use D=3 for screening; verify the genuine transition for the center-blind Higgs reps used here.

4. **Goldstone / spectrum on the lattice.** Confirm the dim_G + (global U(1)) zero-mode structure
   and the massive scalar multiplets dynamically (correlator masses), beyond the tree-level
   Hessian eigenvalues.

5. **Stability of the locking minimum under fluctuations.** Verify that the BFB + stable-locking
   couplings (section 2.5) actually produce the intended H vacuum in a full simulation, including
   the higher-multiplicity singlet directions (e.g. 2T j=6 m=2; the (4,4) flags).

6. **Sigma216 / Sigma648 scalar spectra.** The d=35 and d=64 reps were validated structure-only;
   their full numeric scalar spectra (and m_H) were not computed here.
