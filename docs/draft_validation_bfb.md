# SU(3) bounded-from-below (BFB) conditions for the multi-invariant Higgs quartic

Answers the draft's open `\yy` TODO "check bounded from below for SU(3)" for the
multi-channel scalar potential used in the SU(N) -> discrete-nonabelian-H program.

Code: `src/bfb_su3.cpp` (driver + `check_bfb`), test `test/test_bfb.cpp`.
Build/run:
```
g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -o build/bfb_su3 src/bfb_su3.cpp && OMP_NUM_THREADS=8 ./build/bfb_su3
g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -o build/test_bfb test/test_bfb.cpp && OMP_NUM_THREADS=8 ./build/test_bfb   # 11 passed, 0 failed
```

## Setup

The quartic potential is the channel-projected multi-invariant of `scalar_invariants.hpp`:
```
V_4(Phi) = sum_{rho'} f_{rho'} V_{rho'}(Phi),   V_{rho'}(Phi) = || P_{rho'} M ||_F^2  >= 0,
```
with `M = Phi Phi^dag` (the d x d bilinear) and `P_{rho'}` the projector onto the rho'-th
irreducible channel of `R (x) Rbar`. Those channels are the eigenspaces of the adjoint
quadratic-Casimir superoperator `Chat(M) = sum_a [T^a_R,[T^a_R,M]]`, labelled by their
`C2` eigenvalue -- exactly the draft's `rho'` channels (SU(2): `J`). `V_{rho'} >= 0` always,
being a sum of squared magnitudes `sum_n |Phi^dag Q_{rho'n} Phi|^2`.

For the Sigma(108) Higgs irrep `(2,2) = {4,2}` (d = 27) the channels are
```
C2(rho') = 0, 3, 6, 8, 12, 15, 18, 20, 24      (9 channels)
```
(matches `docs/locking_couplings.md`; `CasimirChannels<3>` returns `0` as `-4.9e-15`).

## Why SU(3) BFB is non-trivial (KEY POINT)

* **SU(2) integer-spin (real) reps.** There the independent quartic invariants reduce to
  **quadratics in the field bilinears**: with a real rep the only invariants are powers of
  `phi^dag phi` and the few channel norms, which for the small integer-`j` cases collapse to
  quadratic forms in a handful of real bilinears. BFB is then a plain quadratic-form
  positivity -- trivial.

* **SU(3).** `V_4 = sum_c f_c ||P_c M||^2` is **genuinely quartic in `Phi`**. It is a
  quadratic form in the entries of `M`, but `M` is *constrained* to the rank-1 positive cone
  `M = Phi Phi^dag`. Positivity of `V_4` over that cone is **copositivity** of the coupling
  set `{f_c}` w.r.t. the channel structure -- NOT plain matrix positivity.

* **Sufficient vs necessary.** `all f_c >= 0` is **SUFFICIENT** (each `V_c >= 0` term-by-term)
  but **NOT necessary**: physical vacuum-alignment minima generally need some `f_c < 0`
  (the Sigma(108) locking set has none negative, but the stable *cone* extends to negative
  `f_c`). How negative an `f_c` may go is bounded by the cross-channel structure
  (`sum_c P_c = 1`, `sum_c C2(c) P_c = Chat`, etc.).

## The BFB criterion and the numeric checker

`V_4` is homogeneous of degree 4: `V_4(t Phi) = t^4 V_4(Phi)`. Hence `V_4` is bounded below
over all of `C^d` **iff** its minimum on the unit sphere `|Phi| = 1` is `>= 0` -- a single
unit vector with `V_4 < 0` already gives the unbounded ray `t Phi, t -> infinity`. So

```
BFB  <=>  min_{|Phi| = 1} V_4(Phi)  >=  0          (copositivity verdict)
```

`check_bfb` minimizes `V_4` on the sphere by many random-restart projected-gradient
descents (Riemannian gradient = free gradient minus its radial part; backtracking line
search; OpenMP over the independent restarts) and declares BFB iff `min >= -eps`
(`eps = 1e-7`). The `argmin` of a non-BFB case IS an unbounded direction.

**Implementation note (speed).** The combined weighted superoperator
`A_f = sum_c f_c P_c` is a single polynomial in `Chat`: since each `P_c = L_c(Chat)`
(Lagrange interpolant), `A_f = poly(Chat)` with `poly(C2(c)) = f_c`. We solve that small
Vandermonde once per `{f}` and apply `A_f M = sum_k coef_k Chat^k M` from the shared Krylov
powers, then `V_4 = <M, A_f M>` and `dV_4/dphibar = 2 (A_f M) Phi`. This replaces
`nc` independent Lagrange chains (`nc*(nc-1)` Casimir applies) by `nc-1` Casimir applies and
matches the direct per-channel sum to `~1e-11` relative (checked in `verify_super`).

## Results

Driver `./build/bfb_su3` on `(2,2) = {4,2}` (Sigma(108) channels):

| coupling set `{f_c}`                       | `min_{|Phi|=1} V_4` | BFB |
|--------------------------------------------|--------------------:|:---:|
| all `f_c = 1`                              |        **+1.000000** | YES |
| `f_last = -1`, rest `0.05`                 |        **-0.490339** | NO  |
| Sigma(108) locking `f` (docs)              |        **+0.052109** | YES |
| `find_stable_couplings` Sigma(108) `f`     |        **+0.052105** | YES |

* **all `f_c = 1` -> `V_4 == 1` exactly.** Because `sum_c P_c = 1` (identity superoperator),
  `V_4 = sum_c ||P_c M||^2 = ||M||_F^2 = (Phi^dag Phi)^2 = 1` for every unit vector. This is
  an exact analytic anchor for the checker (asserted to `1e-6` in the test).
* **`f_last = -1` -> unbounded.** A strongly negative coupling on the highest channel makes
  `V_4 < 0` along that channel's eigendirection; `check_bfb` finds `min ~ -0.49 < 0`.
* **Sigma(108) locking couplings -> BFB.** The couplings that lock SU(3) -> Sigma(108) (a
  stable, positive-semidefinite minimum; `docs/locking_couplings.md`) are a fortiori bounded
  below; `min ~ +0.052 > 0`. Cross-checked with the independently-searched
  `find_stable_couplings` set (same `min`).

### BFB region slice

Hold `f = ` the Sigma(108) locking set and vary `f[0]` (the singlet channel `C2 = 0`) and
`f[last]` (the highest channel `C2 = 24`); `.` = bounded, `X` = unbounded:

```
        y\x  -0.30  -0.10  +0.00  +0.10  +0.30  +0.60
      -0.60   X       X       X       X       X       X
      -0.30   X       X       X       X       X       X
      -0.10   X       X       X       X       X       .
      +0.00   .       .       .       .       .       .
      +0.10   .       .       .       .       .       .
      +0.30   .       .       .       .       .       .
```
(`x = f[0]`, `y = f[last]`.)

**Reading the map (the copositivity structure).** The boundary is set almost entirely by
`f[last]` (the highest-`C2` channel): the region stays bounded down to `f[last] ~ 0`, and
unbounded once `f[last]` goes appreciably negative -- this channel is "load-bearing" for
boundedness. By contrast `f[0]` (the singlet channel) can be pushed strongly negative
(to `-0.3` and beyond) while `V_4` stays bounded, *provided* `f[last] >= 0`. The asymmetry
is the whole content of "some `f_c < 0` are allowed, bounded by cross-channel structure":
boundedness is governed by the channels with the largest Casimir weight, not by the singlet.
(The single `X` at `f[last] = -0.1, f[0] = +0.6` flipping to `.` reflects that a large
positive `f[0]` partially compensates a mildly negative `f[last]` -- copositivity, not
per-channel positivity.)

## Test verdicts (`test/test_bfb.cpp`)

```
(a) all f=1:            vmin = +1.000000e+00  bfb=1      -> BFB; V_4 == ||M||_F^2 == 1
(b) f_last=-1:          vmin = -4.903395e-01  bfb=0      -> NOT BFB (unbounded direction)
(c) Sigma(108) locking: vmin = +5.210864e-02  bfb=1      -> BFB
(c') find_stable_couplings: stable=1  vmin = +5.210522e-02  bfb=1
[test_bfb] 11 passed, 0 failed
```

Asserts: (a) all-positive `f` -> BFB and `V_4 = 1` exactly; (b) negative-dominant `f` ->
NOT BFB; (c) the Sigma(108) locking couplings -> BFB; (c') the independently-searched
stable couplings -> BFB; plus the `(2,2)` channel list matches `0,3,6,8,12,15,18,20,24`.

## Draft cross-check

| draft item | our result | match |
|---|---|---|
| "check bounded from below for SU(3)" (open `\yy` TODO) | BFB = copositivity of `{f_{rho'}}`; checker + region map provided | RESOLVED |
| `(2,2)` `R x Rbar` channels (rho') | `C2 = 0,3,6,8,12,15,18,20,24` (9) | confirmed |
| SU(2) integer-`j` BFB trivial; SU(3) genuinely quartic | confirmed (homogeneous deg-4 copositivity) | confirmed |
| all `f_{rho'} >= 0` sufficient, not necessary | confirmed (region extends to some `f_c < 0`) | confirmed |
| Sigma(108) locking couplings are BFB | `min V_4 = +0.0521 > 0` | confirmed |

No numeric disagreements with the draft were found; the draft left the SU(3) BFB as an
open TODO and this note supplies the criterion (copositivity, equivalently
`min_{|Phi|=1} V_4 >= 0`), a verified numeric checker, and a region map.
