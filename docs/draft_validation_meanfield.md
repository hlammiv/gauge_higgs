# Mean-field strong-weak coupling transition: draft validation + SU(3) extension

Reproduces the tree-level mean-field strong-weak coupling transition (lattice notes
"Mean field method" + "tree level strong-weak coupling transition") for SU(N) gauge
theory frozen onto a finite subgroup H (or kept continuous), and extends the critical
table to SU(3) and the crystal-like Sigma subgroups. Pure group theory + 1-D
root-finding; **no HMC, no lattice simulation**.

- Code: `src/meanfield.cpp` (driver) + `test/test_meanfield.cpp` (asserts).
- Build / run:
  ```
  g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -o build/meanfield src/meanfield.cpp && ./build/meanfield
  g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -Itest -o build/test_meanfield test/test_meanfield.cpp && ./build/test_meanfield
  ```
- Test status: **`[test_meanfield] 48 passed, 0 failed`**.

## Method

Single-link generating log
```
w(alpha) = ln < exp( (alpha/N) Re Tr_fund U ) >,
```
averaged Haar over the continuous group, **uniformly** over the |H| elements for a
discrete subgroup H (Re Tr of the N x N fundamental matrices, from
`finite_subgroups.hpp` `close_group`). Both cases reduce to a finite weighted sum of
nodes `(weight_k, ReTr_k)`, so `w`, `w' = v_B`, ... are cheap log-sum-exp evaluations
(`struct Ensemble` in `meanfield.cpp`).

Closed forms used as cross-checks (both verified in the test):
```
w_SU(2)(a) = ln( 2 I_1(a) / a )                          (I_1 modified Bessel)
w_BT(a)    = ln( 3 + 8 cosh(a/2) + cosh(a) ) - ln 12
```
The 2T (BT) `Re Tr_fund` spectrum is `{-2 (x1), -1 (x8), 0 (x6), +1 (x8), +2 (x1)}`,
so `sum_h exp((a/2) ReTr) / 24 = (3 + 8 cosh(a/2) + cosh(a))/12` -- the group sum and
the closed form agree to 1e-12.

**Continuous SU(3)** uses the maximal-torus diag(e^{i t1}, e^{i t2}, e^{i t3}),
t3 = -t1-t2, with Weyl measure `prod_{i<j} 4 sin^2((t_i - t_j)/2)` and
`Re Tr U = sum_i cos t_i`. Periodic trapezoidal quadrature converges exponentially:
identical to 8 digits already at a 120x120 grid (production uses 360x360). Sanity:
`w_SU(3)(a) -> a^2/36` as a->0 (matches `<(ReTr)^2>_Haar = 1/2`).

### Saddle point (TEMPORAL gauge)
```
N_P/N_l = (D-2)/2                                  (plaquettes per link)
E_P(v)  = 1 - v^4 + (1-v^2) * 2/(D-2)              (mean-field plaquette energy)
v_B     = w'(alpha)                                (mean-link saddle)
alpha   = -beta (N_P/N_l) E_P'(v_B)                (stationarity)
s(v_B)  = w(alpha) - alpha w'(alpha)               (Legendre transform of w)
F_P     = E_P(v_B) - (1/beta)(N_l/N_P) s(v_B)      (free energy per plaquette)
```
The trivial branch is `alpha = 0 -> v_B = 0`, with `F_P^triv = E_P(0) = 1 + 2/(D-2)`.
`beta_c` is the first-order point where the nontrivial branch's `F_P` equals
`F_P^triv` (found by scanning alpha and bisecting the `F_P - F_P^triv` crossing).

## SU(2)-sector criticals -- DRAFT REPRODUCED (temporal gauge)

Computed vs draft reference; all within the 1e-3 target.

### D = 3   (N_P/N_l = 0.500, F_P^triv = 3.0000)
| group       | beta_c  (ref)     | alpha_B  (ref)    | v_B  (ref)        |
|-------------|-------------------|-------------------|-------------------|
| SU(2) cont  | 1.95690 (1.9569)  | 2.43397 (2.4340)  | 0.49822 (0.4982)  |
| 2T (BT)     | 1.93941 (1.9394)  | 3.48834 (3.4883)  | 0.63874 (0.6387)  |
| 2O (BO)     | 1.95665 (1.9566)  | 2.45286 (2.4528)  | 0.50103 (0.5010)  |
| 2I (BI)     | 1.95690 (1.9569)  | 2.43398 (2.4339)  | 0.49822 (0.4982)  |

### D = 4   (N_P/N_l = 1.000, F_P^triv = 2.0000)
| group       | beta_c  (ref)     | alpha_B  (ref)    | v_B  (ref)        |
|-------------|-------------------|-------------------|-------------------|
| SU(2) cont  | 1.68173 (1.6817)  | 4.71901 (4.7190)  | 0.70430 (0.7043)  |
| 2T (BT)     | 1.53740 (1.5374)  | 7.73094 (7.7309)  | 0.92604 (0.9260)  |
| 2O (BO)     | 1.67377 (1.6738)  | 5.05231 (5.0523)  | 0.73029 (0.7303)  |
| 2I (BI)     | 1.68171 (1.6817)  | 4.72008 (4.7201)  | 0.70438 (0.7044)  |

All 8 numbers reproduced to <= 1.1e-4 of the draft table (well inside 1e-3).
2I (|H| = 120) is so dense in SU(2) that its critical is numerically degenerate with
the continuous SU(2) value -- its discrete frozenness only shows up at much larger
alpha (see below).

## SU(3)-sector criticals -- NEW (temporal gauge)

N = 3, exponent `alpha/3`, `Re Tr` of the 3x3 group elements (`Sigma` via
`close_group`); continuous SU(3) via the maximal-torus Weyl quadrature.

### D = 3   (N_P/N_l = 0.500, F_P^triv = 3.0000)
| group        |  beta_c  |  alpha_B  |   v_B    |   1 - v_B   |
|--------------|----------|-----------|----------|-------------|
| SU(3) cont   | 6.15144  | 12.59143  | 0.69203  | 3.080e-01   |
| Sigma(108)   | 3.11794  | 12.37905  | 0.99627  | 3.727e-03   |
| Sigma(216)   | 3.58079  | 14.24242  | 0.99718  | 2.825e-03   |
| Sigma(648)   | 4.30518  | 16.95898  | 0.99236  | 7.644e-03   |
| Sigma(1080)  | 4.64650  | 18.31217  | 0.99259  | 7.407e-03   |

### D = 4   (N_P/N_l = 1.000, F_P^triv = 2.0000)
| group        |  beta_c  |  alpha_B  |   v_B    |   1 - v_B   |
|--------------|----------|-----------|----------|-------------|
| SU(3) cont   | 5.07223  | 17.11209  | 0.77083  | 2.292e-01   |
| Sigma(108)   | 2.34019  | 14.00158  | 0.99879  | 1.210e-03   |
| Sigma(216)   | 2.68704  | 16.09162  | 0.99919  | 8.143e-04   |
| Sigma(648)   | 3.23397  | 19.28219  | 0.99731  | 2.693e-03   |
| Sigma(1080)  | 3.48993  | 20.82708  | 0.99769  | 2.308e-03   |

(Group orders confirmed by `close_group`: |108|, |216|, |648|, |1080|.)

Trends: discrete Sigma groups transition at noticeably **smaller** beta_c than
continuous SU(3) (frozen links are energetically cheaper to align), and their saddle
mean-link sits almost on the unit circle (`v_B ~ 0.99`), i.e. the link has essentially
locked onto a group element -- whereas continuous SU(3) only reaches `v_B ~ 0.69`
(D=3) / `0.77` (D=4). beta_c grows monotonically with |H| across the Sigma series
(108 < 216 < 648 < 1080), as the larger groups approximate SU(3) more closely.

## Discrete frozenness signature

`1 - v_B(alpha) = 1 - w'(alpha)` as alpha grows: **continuous G decays as a power law,
discrete H decays exponentially** (the link freezes onto the nearest group element,
controlled by the gap between the largest `Re Tr` = N and the next distinct value).

| alpha | SU(2)    | 2T (BT)  | 2I (BI)  | SU(3)    | Sigma(108) |
|-------|----------|----------|----------|----------|------------|
|  4.0  | 3.42e-01 | 3.08e-01 | 3.42e-01 | 7.20e-01 | 6.10e-01   |
|  8.0  | 1.81e-01 | 6.56e-02 | 1.80e-01 | 4.64e-01 | 7.95e-02   |
| 12.0  | 1.22e-01 | 9.76e-03 | 1.14e-01 | 3.22e-01 | 4.86e-03   |
| 16.0  | 9.22e-02 | 1.34e-03 | 7.09e-02 | 2.45e-01 | 3.07e-04   |
| 20.0  | 7.40e-02 | 1.82e-04 | 4.01e-02 | 1.97e-01 | 2.04e-05   |
| 30.0  | 4.96e-02 | 1.22e-06 | 7.17e-03 | 1.32e-01 | 2.49e-08   |

- SU(2): `1 - v_B ~ 3/(2 alpha)` -- at alpha=30, 4.96e-2 ~ 0.05 (power law).
- SU(3): also power law (1.3e-1 at alpha=30).
- 2T (BT): `1 - v_B ~ e^{-alpha/2}` (gap = N - 1 = 1 between ReTr = 2 and ReTr = 1,
  times alpha/N = alpha/2) -- drops ~150x over alpha = 20 -> 30 (exponential).
- Sigma(108): even steeper exponential (2.5e-8 at alpha=30).
- 2I (BI) is intermediate: being the densest SU(2) subgroup (|H| = 120) it tracks the
  continuous SU(2) power law until large alpha, where its own discrete gap finally
  drives an exponential collapse -- consistent with its near-degenerate critical.

This exponential frozenness is the mean-field face of the deconfinement / link-locking
that distinguishes a discrete gauge group from its continuous parent.

## Test coverage (`test/test_meanfield.cpp`, 48 checks)
- |2T| = 24, |2O| = 48, |2I| = 120.
- `w_BT` closed form == group sum (6 alpha values, 1e-12).
- `w_SU(2)` Bessel form == torus quadrature (5 alpha values, 1e-6).
- `w(0) = 0`, `w'(0) = v_B(0) = 0` (symmetric ensembles).
- All 8 SU(2)-sector criticals (beta_c, alpha_B, v_B) vs the draft table at 1e-3.

## Findings / notes
- No discrepancies with the draft. Every SU(2)-sector entry reproduces to <= 1.1e-4,
  comfortably inside the stated 1e-3 tolerance.
- The SU(3) sector numbers above are **new** (not in the draft). They are produced by
  the identical saddle-point machinery validated against the SU(2) table, with the
  continuous-SU(3) `w` cross-checked by exponential quadrature convergence and the
  small-alpha `a^2/36` limit, and the Sigma `w` being exact finite group sums.
