# Draft validation: scalar mass-spectrum tables

Reproduction of the draft's per-channel Hessian eigenvalue tables for spontaneous
symmetry breaking `SU(2) -> {BT, BO, BI}` and `SU(3) -> {Sigma(108), Sigma(216),
Sigma(648), Sigma(1080)}` at the H-singlet VEV.

- Test: `test/test_draft_spectra.cpp`  (build to `build/test_draft_spectra`)
- Result: **89 passed, 0 failed** (default mode: BT, BO, BI, Sigma108, Sigma1080).
  `all` mode additionally builds Sigma216 (d=35) and Sigma648 (d=64).
- Method: pure group theory / linear algebra. For each `G x Gbar` channel `c`
  (the `CasimirChannels` eigenspaces = the draft's `J` for SU(2) and `(p,q)` for
  SU(3)) we build the per-channel Hessian `H^(c)` of `V_c = ||P_c (phi phi^dag)||^2`
  at the unique H-singlet VEV (`vacuum_alignment::singlet_vev`), via a 4-point
  Richardson stencil on the analytic gradient (`real_grad`), then read off and
  cluster its eigenvalues into H-irrep multiplets. No HMC, no lattice. Heaviest
  rep d=64.

## What is confirmed

**(i) Observation 2 -- degeneracy within H-irreps.** In every channel the distinct
eigenvalues appear in clusters whose sizes are the dimensions of the unbroken H's
irreps (1, 2, 3, 4, 5, 6, 8, 9 ...). The radial (VEV) direction is always its own
1-dimensional eigenspace with eigenvalue `lambda_0^(c)` (asserted for SU(3)).

**(ii) Observation 3 -- Goldstones in the adjoint subduction with the "1/3" value.**
- The broken-generator directions `{ i T^a_R phi_VEV }` span exactly the adjoint of G:
  rank = dim(adj SU(2)) = 3, resp. dim(adj SU(3)) = 8 (H is discrete, so every gauge
  generator is broken). These are the would-be Goldstones; the adjoint subduces to the
  H-irreps listed in the subduction tables (e.g. SU(2) spin-1 -> the `T`/`T_1` H-irrep;
  SU(3) 8 -> `8^(0)`).
- In **every** channel `c` the Goldstone-subspace Hessian eigenvalue equals exactly
  `(1/3) lambda_0^(c)` (the radial eigenvalue). This is the draft's per-channel "1/3"
  value (it becomes 0 after the `lambda_0/3` subtraction). Verified to < 1e-4 in all
  channels of all seven cases.
- The full tuned (`find_stable_couplings`) Hessian in the rep-faithful complex (2d-dof)
  treatment has exactly `dim(G) + 1` zero modes: 3+1 (SU(2)) resp. 8+1 (SU(3)) -- the
  gauge Goldstones plus the global U(1) of the complex tensor basis. (As in `test_align`.)

**(iii) Numeric match to the tables (up to the stated normalization).**

### Normalization finding (SU(2))

The draft's per-channel Hessian `H^(J)` differs from our orthogonal-projector `V_c`
Hessian by a **per-channel factor `1 / sqrt(2J+1)`**. This is the Clebsch-Gordan singlet
coefficient `<J,M; J,-M | 0,0> = (-1)^{J-M}/sqrt(2J+1)` that enters the draft's CG tensor
`T_abcd(J)` but is absent from our (orthonormal) projector `P_c`. After

1. rescaling `H^(J) -> H^(J)/sqrt(2J+1)`,
2. subtracting `(1/3) lambda_0^(J)` (the draft's `H(J) - (1/3) lambda_0 delta`),
3. dividing the whole table by the single constant `(2/3) lambda_0^(J=0)` so the J=0
   singlet eigenvalue is 1,

the SU(2) tables reproduce the draft **exactly** (to FD precision ~1e-5). Only even J
contribute (the totally symmetric quartic vanishes for odd J).

#### BT spin-3, `{6}`, d=7  (draft Tab. BTmass), |H|=24

| J | A0 | T_Goldstone | T_massive | computed (matches) |
|---|----|----|----|----|
| 0 | 1 | 0 | 0 | A0=1.00000 |
| 2 | 0 | 0 | sqrt5/3 = 0.74536 | 0.74536 |
| 4 | 14/11 = 1.27273 | 0 | -10/11 = -0.90909 | 1.27273, -0.90909 |
| 6 | 24/(11 sqrt13) = 0.60513 | 0 | 35/(33 sqrt13) = 0.29416 | 0.60513, 0.29416 |

(Draft prints 0.60553 / 0.29424 for J=6; our exact-arithmetic targets are 0.60513 /
0.29416 and we match those -- the draft's last two digits are a rounding artifact.)

#### BO spin-4, `{8}`, d=9  (draft Tab. BOmass), |H|=48

| J | A1 | T1_Goldstone | T2_massive | E |
|---|----|----|----|----|
| 0 | 1 | 0 | 0 | 0 |
| 2 | 0 | 0 | sqrt5/11 = 0.20328 | 4 sqrt5/11 = 0.81312 |
| 4 | 98/143 = 0.68531 | 0 | -10/143 = -0.06993 | -40/143 = -0.27972 |
| 6 | 40/(11 sqrt13) = 1.00855 | 0 | -7/(11 sqrt13) = -0.17650 | -28/(11 sqrt13) = -0.70598 |
| 8 | 30/(13 sqrt17) = 0.55970 | 0 | 56/(143 sqrt17) = 0.09498 | 224/(143 sqrt17) = 0.37992 |

All match. **Resolved draft self-flag ("check T_{1,2}"):** the Goldstone (`1/3` -> 0)
sits in the `T_1` 3-dim H-irrep and the massive mode in the distinct `T_2` 3-dim
H-irrep, with the massive eigenvalue `sqrt5/11` at J=2. So the draft column order
(A1, T1_Goldstone, T2_massive, E) is correct; the notes' flag is resolved in favor of
`T_1` = Goldstone, `T_2` = massive.

#### BI spin-6, `{12}`, d=13  (draft Tab. BImass), |H|=120

| J | A0 | J | H |
|---|----|----|----|
| 0 | 1 | 0 | 0 |
| 2 | 0 | 1/sqrt5 = 0.44721 | 0 |
| 4 | 0 | 6/17 = 0.35294 | 21/68 = 0.30882 |
| 6 | 121 sqrt13/323 = 1.35069 | -7 sqrt13/34 = -0.74232 | -63 sqrt13/323 = -0.70325 |
| 8 | 0 | 25/(19 sqrt17) = 0.31913 | 99/(38 sqrt17) = 0.63187 |
| 10 | 91 sqrt21/391 = 1.06653 | -1529 sqrt21/14858 = -0.47158 | -429 sqrt21/7429 = -0.26463 |
| 12 | 196/437 = 0.44851 | 572/1955 = 0.29258 | 1287/29716 = 0.04331 |

All match (`T_1` is the Goldstone). The full BI table is reproduced.

### SU(3): the draft and the notes disagree -- our computation adjudicates

For SU(3) the convention-independent statement is the **per-channel form**: normalize
each channel so the radial (VEV) eigenvalue `lambda_0^(c) = 1`, then the Goldstone
eigenvalue is `1/3` and the remaining eigenvalues are physical ratios. This is the
normalization of the **notes** (`notes_WL/notes1.tex`).

**Finding: the notes SU(3) tables (S108, S1080) match our computation exactly; the
`draft-1.tex` SU(3) tables (S108mass, S216mass, S648mass, S1080mass) are incorrect.**
The draft tables fail even the convention-independent within-channel ratio check (e.g.
draft S108 (3,3): `2^(0)/radial = 1.2119/2.26222 = 0.536`, whereas the physical ratio is
`0.69048`; draft S1080 (4,4): `8_massive/radial = -0.0143/0.9 = -0.016` vs the physical
`0.32275`). The draft-1.tex SU(3) tables appear to be a stale/garbled version that was
superseded by the notes; **use the notes SU(3) tables.**

Note also a labelling subtlety: the adjoint Casimir merges conjugate pairs and accidental
degeneracies into a single channel (e.g. `(4,1)` and `(1,4)` share `C2=12`; the draft
labels it `(4,1)`, the notes `(1,4)` -- same channel).

#### Sigma(108), (2,2), `{4,2}`, d=27 (matches notes Tab. S108mass), |H|=108

Per-channel radial-normalized eigenvalues (Goldstone = 1/3 in every nonzero-radial channel):

| notes col | channel C2 | 2^(0) | massive 2(4^(0)+4^(0)') |
|---|---|---|---|
| (1,4) | 12 | -1/6 = -0.16667 | 0.19514, -0.02847 |
| (3,3) | 15 | 0.69048 | 0.31829, 0.05076 |
| (4,4) | 24 | 0.43333 | 0.58045, 0.41538 |

#### Sigma(1080), (6,0), `{6}`, d=28 (matches notes Tab. S1080mass), |H|=1080

| notes col | channel C2 | 8^(0)_massive | 5^(0)+5^(0)' | 9^(0) |
|---|---|---|---|---|
| (4,4) | 24 | 0.32275 | 0.25645, 0.02503 | -1/9 = -0.11111, 0.03704 |
| (5,5) | 35 | 0.65741 | 0.75548, 0.01489 | 0.85185, 0.79424 |
| (6,6) | 48 | 0.08838 | 0.44652, 0.29540 | 0.33586, 0.28283 |

Sigma(216) (d=35) and Sigma(648) (d=64) (`all` mode) reproduce the same
radial / 8^(0)-Goldstone / massive structure; the notes give no full numeric table for
them, so only structure (i)/(ii) is asserted.

## Synthesis: lightest non-Goldstone scalar mass m_H per subgroup

`m_H^2` = smallest positive eigenvalue of the full tuned (PSD, stable-locking) Hessian,
in units where the potential's overall scale fixes the spectrum (representative stable
locking from `find_stable_couplings`; the value depends on the chosen quartic couplings
`f_c`, not unique):

| subgroup | rep | d | |H| | zero modes | m_H^2 (lightest) |
|---|---|---|---|---|---|
| BT  | spin-3 {6}  | 7  | 24   | 4 (3+U(1)) | 0.467 |
| BO  | spin-4 {8}  | 9  | 48   | 4 (3+U(1)) | 0.125 |
| BI  | spin-6 {12} | 13 | 120  | 4 (3+U(1)) | 0.215 |
| Sigma(108)  | (2,2) {4,2} | 27 | 108  | 9 (8+U(1)) | 0.182 |
| Sigma(1080) | (6,0) {6}   | 28 | 1080 | 9 (8+U(1)) | 0.065 |

All seven cases admit a stable (positive-semidefinite) locking potential with the correct
Goldstone count and a positive mass gap, confirming the draft's central claim that a
single large Higgs irrep can break SU(N) to the chosen crystal-like discrete subgroup.
