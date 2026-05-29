# Draft validation: gauge-boson mass matrix M_ab

Validates the draft sections "Gauge bosons spectra" and "Enhanced symmetries" (eq. M_ab).

**Test:** `test/test_draft_gaugemass.cpp` -> `build/test_draft_gaugemass`
**Build/run:**
```
g++ -std=c++20 -O3 -march=native -fopenmp -Isrc -Itest -o build/test_draft_gaugemass test/test_draft_gaugemass.cpp
./build/test_draft_gaugemass
```
**Result:** `[test_draft_gaugemass] 60 passed, 0 failed`. No HMC / lattice; pure group theory + linear algebra on rep dim <= 64.

## What is computed

With `phi_B` the H-singlet VEV (image of the group-averaging projector `(1/|H|) sum_h D^(sigma)(h)`, from `vacuum_alignment::singlet_vev`) and `T^(sigma)_a` the irrep generators (from `GeneralRep<N>`):

```
M_ab = Re( phi_B^dagger T^(sigma)_b T^(sigma)_a phi_B ).
```

Because each `T^a` is Hermitian, `M_ab = Re( (T_b phi)^dagger (T_a phi) )`, which is real-symmetric. The H-singlet VEV is unit-normalized (`|phi_B| = 1`).

## Claims verified (all 7 subgroups, tight tolerances)

**(i) M is proportional to the identity** — off-diagonal ~ 0, all diagonal entries equal, all eigenvalues equal and strictly positive. Hence M is rank-full over all `dim(G) = N^2-1` generators: every gauge boson is degenerate and the **continuous gauge symmetry is completely broken** (no surviving continuous subgroup; only the discrete H remains).

| subgroup | max off-diag | diag spread | eig spread |
|---|---|---|---|
| BT (2T) | 0 | 0 | 0 |
| BO (2O) | 5.4e-35 | 2.7e-15 | 2.7e-15 |
| BI (2I) | 1.6e-14 | 2.0e-13 | 2.2e-13 |
| Sigma108 | 3.3e-16 | 1.6e-15 | 1.6e-15 |
| Sigma216 | 1.7e-16 | 4.9e-15 | 4.9e-15 |
| Sigma648 | 1.9e-15 | 2.1e-14 | 2.1e-14 |
| Sigma1080 | 1.2e-15 | 8.0e-15 | 8.0e-15 |

**(ii) Tr M = C2(sigma)** — the quadratic Casimir of the Higgs irrep. This is an identity:
`Tr M = sum_a phi^dag (T^a)^2 phi = phi^dag (sum_a (T^a)^2) phi = C2(sigma) |phi|^2 = C2(sigma)`.
Matched to 1e-7 for all 7.

**(iii) Common eigenvalue / gauge mass** — `M = c I` with `c = Tr M / dim_G = C2(sigma)/dim_G`. The physical gauge mass squared (kappa g^2 = 1 units) carries the draft's factor of 2:
```
m_A^2 = (2/dim_G) kappa g^2 C2(sigma)  = 2 * c.
```
For SU(3) this is `m_A^2 = (1/12)(p^2+q^2+pq+3p+3q)` (dim_G = 8); for SU(2) `m_A^2 prop to C2(j)=j(j+1)` (dim_G = 3).

## Conventions

- Generators normalized `Tr(T^a T^b) = (1/2) delta^{ab}` (`algebra.hpp`). In this normalization the SU(3) quadratic Casimir is `C2(p,q) = (1/3)(p^2+q^2+pq+3p+3q)` and the SU(2) spin-j Casimir is `C2(j) = j(j+1)`.
- SU(2): GeneralRep<2>({2j}); SU(3): (p,q) = Young rows [p+q, q]. So spin-3={6}, spin-4={8}, spin-6={12}; (2,2)={4,2}, (4,1)={5,1}, (3,3)={6,3}, (6,0)={6}.
- **Factor of 2 (NOT a discrepancy):** the *matrix* `M_ab` defined above has common diagonal `C2/dim_G`; the *physical* gauge mass^2 is `m_A^2 = (2/dim_G) C2 = 2 * (common diagonal)`. The 2 is the standard kinetic-normalization factor in the draft's `m_A^2 = (2/dim_G) kappa g^2 C2`. Both forms are reported below.

## SU(3) draft bracket cross-check

The draft brackets `p^2+q^2+pq+3p+3q` are reproduced exactly:
`(2,2)->24`, `(4,1)->36`, `(3,3)->45`, `(6,0)->54`, giving `m_A^2 = (1/12)*bracket = 2, 3, 3.75, 4.5`.

## Synthesis table (kappa g^2 = 1) — gauge mass m_A for type-I/II comparison vs scalar m_H

| subgroup | N | dim_G | dim_R | \|H\| | C2(sigma) | common eig (C2/dim_G) | **m_A^2 = 2*eig** | **m_A** |
|---|---|---|---|---|---|---|---|---|
| BT (2T, spin-3)   | 2 | 3 | 7  | 24   | 12 | 4.00000   | 8.00000  | 2.82843 |
| BO (2O, spin-4)   | 2 | 3 | 9  | 48   | 20 | 6.66667   | 13.33333 | 3.65148 |
| BI (2I, spin-6)   | 2 | 3 | 13 | 120  | 42 | 14.00000  | 28.00000 | 5.29150 |
| Sigma108 (2,2)    | 3 | 8 | 27 | 108  | 8  | 1.00000   | 2.00000  | 1.41421 |
| Sigma216 (4,1)    | 3 | 8 | 35 | 216  | 12 | 1.50000   | 3.00000  | 1.73205 |
| Sigma648 (3,3)    | 3 | 8 | 64 | 648  | 15 | 1.87500   | 3.75000  | 1.93649 |
| Sigma1080 (6,0)   | 3 | 8 | 28 | 1080 | 18 | 2.25000   | 4.50000  | 2.12132 |

m_A is in units where kappa g^2 = 1 (and the H-singlet VEV is unit-normalized; for a VEV of length v multiply M, hence m_A^2, by v^2).

## Findings

No discrepancies. All draft claims (i)/(ii)/(iii) and the SU(3) brackets reproduce exactly. The only convention point worth flagging for the synthesis: the draft's `m_A^2 = (2/dim_G) kappa g^2 C2` carries an explicit factor of 2 relative to the raw mass-matrix diagonal `M_aa = C2/dim_G`; the table reports both `common eig` (= M_aa) and `m_A^2 = 2*eig` so the synthesis can use the physical mass directly.
