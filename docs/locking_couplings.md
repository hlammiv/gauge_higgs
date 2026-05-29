# Locking couplings f_c for SU(N) -> discrete nonabelian H

Stable multi-invariant scalar-potential couplings that break SU(N) to a crystal-like
discrete subgroup H, computed by `test/test_align.cpp` (build `build/test_align`,
run with `all` for the heavy SU(3) reps). Potential:
```
V(phi) = -mu^2 (phi^dag phi) + sum_c f_c V_c[phi],   V_c = Tr(M^dag P_c M), M = phi phi^dag,
```
where P_c projects the bilinear onto the c-th adjoint-Casimir channel of R (x) Rbar
(channel label = its C2 eigenvalue). Couplings live on the simplex (sum f_c = 1).

Each entry below is ONE validated stable point (from a random simplex search): the VEV
is the H-singlet (image of the group-averaging projector, H-invariance ~1e-15), the
Hessian is positive semidefinite (no negative modes), with exactly (N^2-1) gauge
Goldstones + 1 global-U(1) Goldstone as zero modes, and the massive scalars fall into
H-irrep multiplets (draft Observation 2). The f_c are NOT unique — any point in the
stable cone locks H; these are concrete starting points for HMC runs. NOT yet checked:
that this stratum is the GLOBAL minimum (vs other subgroups) — a further scan.

| H | rep | d | channels C2 (order of f_c) | mu^2 | mass gap | massive multiplets (m^2 x deg) |
|---|---|---|---|---|---|---|
| 2T (spin-3) | {6} | 7 | 0,2,6,12,20,30,42 | 0.113 | 0.453 | A0(x1)+ T triplets x3 |
| 2O (spin-4) | {8} | 9 | 0,2,6,12,20,30,42,56,72 | 0.108 | 0.126 | x3,x2,x3,x3,x1,x2 (T,E,A) |
| 2I (spin-6) | {12} | 13 | 0,2,6,12,20,30,42,56,72,90,110,132,156 | 0.065 | 0.213 | 2I multiplets (mostly x1; some x2,x3,x4) |
| Sigma(108) | (2,2)={4,2} | 27 | 0,3,6,8,12,15,18,20,24 | 0.104 | 0.161 | octets x8 + x1,x2 |
| Sigma(216) | (4,1)={5,1} | 35 | 0,3,6,8,12,15,20,24,30,35 | 0.109 | 0.143 | x8,x4,x8,x16,x16,x8,x1 |
| Sigma(1080) | (6,0)={6} | 28 | 0,3,8,15,24,35,48 | 0.027 | 0.094 | (9 zero modes; positive massive) |

### f_c values (per channel, in the C2 order above)
- **2T**:        0.1287 0.1548 0.1835 0.2399 0.0056 0.1745 0.1130
- **2O**:        0.0562 0.2045 0.2230 0.2206 0.0444 0.0740 0.0066 0.0340 0.1366
- **2I**:        0.1182 0.0819 0.0353 0.0274 0.1550 0.0828 0.0157 0.0702 0.1083 0.1143 0.0189 0.1124 0.0595
- **Sigma(108)**:  0.2782 0.0930 0.1029 0.1256 0.0071 0.0101 0.0443 0.2148 0.1240
- **Sigma(216)**:  0.0213 0.0726 0.0376 0.2013 0.1806 0.0053 0.1992 0.0513 0.0476 0.1831
- **Sigma(1080)**: 0.1079 0.0906 0.3809 0.3968 0.0060 0.0028 0.0152

6 of 7 crystal-like subgroups done. (2I unblocked by the direct symmetric-rep basis for
single-row Young diagrams; its channel decomposition uses double-reorthogonalized Lanczos.)

### Remaining
- **Sigma(648) ((3,3), d=64)**: alignment works but the per-channel Hessian precompute is
  O(N^2 d^3) heavy at d=64 (~tens of minutes serial). Needs the cached/optimized rep path.
