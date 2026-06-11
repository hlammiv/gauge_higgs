# DEFINITIVE gauge-sector observable protocol — U(1) + charge-q Higgs

Status: settled. This supersedes all prior internal disagreement. Scope: the three phases of the
4D compact U(1) + charge-q fixed-length Higgs model (Confined / Coulomb / Higgs), for q = 1, 2, 4, 6, 8,
including the q≥5 intermediate Coulomb wedge (task 26) and the triple point (task 28).

Convention: our `beta` = gauge coupling = Bowler `beta_p` = Damgaard-Heller `beta` = Fradkin-Shenker `K`.
Our `kappa` = Higgs hopping = Bowler `beta_L` = D-H `kappa` = FS `P` (FS axis is transposed: FS beta = our kappa).
`q` = charge = Bowler/D-H `Q`. Fixed-length |phi|=1 = lambda→inf (no `<rho^2>` length observable; that is a
variable-length / radial-route diagnostic only).

---

## 0. EXECUTIVE SUMMARY (read this; the rest is justification)

THE MINIMAL SUFFICIENT IDENTIFY SET IS THREE OBSERVABLES, NOT TWO:

    { m_gamma (photon mass, structure factor) , sigma_1 (charge-1 string tension) , rho_M (monopole density) }

with `|<P_1>|` (charge-1 Polyakov modulus) as the cheap, single-config cross-check on sigma_1, and `sigma_q`
(charge-q string tension) as the screened negative control. The pair {m_gamma, rho_M} ALONE is INSUFFICIENT for
q≥2 — that is the resolved conflict.

THE PHASE SIGNATURE TRIPLE  (m_gamma, sigma_1, rho_M):

| Phase                              | m_gamma | sigma_1 | rho_M | |<P_1>| |
|------------------------------------|---------|---------|-------|--------|
| Confined                           | massive | area>0  | HIGH  | ~0     |
| Coulomb (incl. q≥5 wedge)          | **0**   | 0/perim | low   | >0     |
| Higgs, deep-Higgs Z_q-confined (q≥2)| massive | **area>0** | low | ~0     |
| Higgs, Z_q-deconfined (large beta) | massive | 0/perim | low   | >0     |

The deep-Higgs Z_q-confined row is the one {m_gamma, rho_M} cannot see (both read the same as ordinary Higgs);
ONLY sigma_1 / |<P_1>| flags it. For q=1 the bottom three rows collapse — sigma_1=0 everywhere in the matter
region (charge-1 IS the condensate charge, screened), so Confined and Higgs are one analytically-connected phase
(Fradkin-Shenker), and there is no triple point of the q≥2 type. q=1 sigma_1 is the built-in negative control.

THE CONFLICT RULING:  `u1-deep-higgs-zq-blindspot.md` is RIGHT. `u1-coulomb-locator-photon-mass.md`'s line
"m_gamma + rho_M separate all three phases" is FALSE for q≥2 and is retired (it is accidentally correct only for
q=1 / outside the deep-Higgs Z_q corner). Reason: m_gamma and rho_M are not independent (M_gamma^2 ∝ n_monopole),
so they collapse toward one axis and jointly separate {massless Coulomb} from {massive, monopole-screened} but
CANNOT split the massive side into Higgs vs Z_q-confined. sigma_1 is mandatory because charge-1 is unscreened for
q≥2 (1 ≠ 0 mod q) and is the only observable that is area-law in the deep-Higgs Z_q corner.

THE WIRING LIST (mechanical port onto src/u1_frozen.cpp; no new physics):
  1. `wilson_grids<kDim>(s.th, s.lat, q, Rmax, w1, wq)` + `CreutzJackknife` jk1/jkq + `reliableR`/`plateau_sigma_jack`
     → sigma_1 ± err and sigma_q ± err.  (port from u1_scan.cpp:52-90, 287-307)
  2. `polyakov_abs<kDim>(s.th, s.lat, 1)` accumulated in a `Stats` → <|P_1|> ± err and chi_{P1}=Vsp·Var.  (u1_scan.cpp:96-109)
  3. `monopole_density<kDim>(s.th, s.lat)` accumulated in a `Stats` → rho_M ± err.  (u1/monopole.hpp:138)
  4. `photon_structure_factor<kDim>(s.th, s.lat, mom)` per config → `photon_mass_fit` → m_gamma ± err.
     ONLY trustworthy at spatial extent L_s ≥ 16 (anisotropic 16^3 × L_t); emit but gate interpretation.  (u1/photon_structure.hpp)
No new observable needs writing: |<P_1>| already exists as `polyakov_abs`. Do NOT use `polyakov(...,m)` (returns
<cos(line)>≈0, not the modulus) and do NOT use the temporal correlator in `photon_mass.hpp` (sum-rule + doubler;
the `m_gamma`/`m_gamma_cosh` CSV columns are the broken, inverted ones — drop them).

---

## 1. THE TABLE (minimal sufficient gauge observables × the three phases)

"area" = positive Creutz/Wilson string tension (confining). "perim" = zero / perimeter law (deconfined).
"massless" = m_gamma = 0. For q≥2 the Confined column means "any confinement"; its q-dependence is split into
the small-kappa monopole-driven part and the large-kappa Z_q-flux-driven part (the deep-Higgs tongue).

| Observable | Confined (small beta) | Coulomb (large beta, small kappa; + q≥5 wedge to kappa=inf) | Higgs (large kappa) | Boundary it identifies |
|---|---|---|---|---|
| **m_gamma** — photon mass via magnetic structure factor R^-1 = a + b·phat^2, m^2 = a/b, p_0=0 (photon_structure.hpp) | massive (>0); m ~ -ln(beta) | **= 0 (MASSLESS) — the UNIQUE massless phase** | massive (>0); m^2 = kappa/beta (Higgs mechanism) | **Coulomb ↔ {Confined, Higgs}** (the only massless/massive line). Flags the q≥5 wedge as the m_gamma=0 strip persisting to large kappa. Needs L_s ≥ 16. |
| **sigma_1** — charge-1 Creutz/Wilson string tension (wilson_grids g1 → CreutzJackknife → plateau_sigma_jack) | area > 0 | 0 / perim | **area > 0 in deep-Higgs Z_q-confined (q≥2, small beta); 0/perim once charge-1 deconfines (large beta).** q=1: 0/perim everywhere (screened control) | **{Confined incl. deep-Higgs Z_q} ↔ {Coulomb, Higgs}.** THE deep-Higgs Z_q line that climbs to kappa=inf (q=2 → Z_2 self-dual beta=0.4407). The ONLY observable not blind there. |
| **rho_M** — monopole density, DeGrand-Toussaint (monopole.hpp) | **HIGH** (monopole plasma) | low (dilute) | low (charge-q condensate Higgs-SCREENS U(1) monopoles → collapses ~0.3 by kappa~0.5, all q) | **CORROBORATOR only.** Sub-classifies confinement (HIGH = monopole-driven, low = Z_q-flux-driven). BLIND to the deep-Higgs Z_q tongue. UV-noisy below beta~1. Never a phase identifier on its own. |
| **\|<P_1>\|** — charge-1 Polyakov modulus, polyakov_abs (u1_scan.cpp:96) | ~0 | > 0 | ~0 in Z_q-confined; > 0 in Z_q-deconfined / large beta. q=1: screened control | Cheap cross-check on sigma_1 (same cut). chi_{P1}=Vsp·Var peak co-locates the Z_q-deconfinement line. Finite-V biased upward; corroborates, does not replace sigma_1. |
| **sigma_q** — charge-q string tension, screened control (wilson_grids gq) | area > 0 | 0 | 0 (screened: q\|q always) | NEGATIVE control. sigma_q→0 while sigma_1>0 certifies the residual is **Z_q, not trivial**. |
| **chi_plaq = Var(A)/V** — gauge specific heat (already in u1_frozen) | peak at confinement/freezing line | — | peak at the q-dependent deep-Higgs Z_q gauge transition | LOCATOR (q-dependent). Locates the gauge/confinement line + Z_q line position. Anchor top at pure-U(1) beta_c≈1.01. |
| **chi_link = Var(B)/V, <cos>_link** — matter specific heat (already in u1_frozen) | — | — | monotone CROSSOVER (no peak; q-BLIND, identical q=2..8) | LOCATOR-ONLY. Marks where matter energy turns over. **Does NOT identify** the Higgs line (Elitzur: no local Higgs order parameter). |
| photon_mass.hpp temporal correlator (`m_gamma`/`m_gamma_cosh` columns) | INVALID | INVALID (flux sum-rule forces negative tail; F_{0i} slice-straddle; p_t=pi doubler) | INVALID | **DO NOT USE.** Reads false-massless inside the Higgs region. Drop the columns. |

q-dependence of the topology (the one fact that sets the whole diagram):

- **q = 1:** TWO phases only. Confined and Higgs are analytically connected (FS; Bowler fig.1: no transition for
  beta≤0.5 or kappa≥1). sigma_1 = 0 throughout the matter region. NO triple point.
- **q = 2, 3, 4:** kappa=inf residual is pure Z_q with ONE transition (Creutz-Jacobs-Rebbi; no Coulomb at the
  edge). FINITE triple point (beta_t, kappa_t). The deep-Higgs Z_q-confined tongue runs to kappa=inf, ending at
  the Z_q self-dual point (q=2: beta=0.4407 = ½ln(1+√2)). sigma_1 carries this line; rho_M and m_gamma are blind to it.
- **q ≥ 5:** kappa=inf residual Z_q has TWO transitions with a "Z_q-symmetric, massless, deconfining"
  intermediate phase between them (CJR; D-H fig.2/4). NO finite triple point: the massless Coulomb WEDGE
  (m_gamma=0 AND sigma_1=perim) persists all the way to kappa=inf, between the lower Z_q transition (beta~1)
  and the upper one (beta ~ O(q^2): ~36 for q=6, ~64 for q=8). The Higgs/ordered phase is pushed into the
  upper-right (beta > O(q^2)) corner and shrinks as q grows.

---

## 2. THE RULING (m_gamma + rho_M  vs  sigma_1)

DEFINITIVE: `u1-deep-higgs-zq-blindspot.md` is CORRECT and load-bearing. The claim in
`u1-coulomb-locator-photon-mass.md` that "m_gamma + rho_M separate all three phases" is FALSE for q≥2 and is
RETIRED. Both memories are right about their own observable's behavior; the photon-mass memory over-generalized
the q=1 topology to all q.

WHY (paper- and code-backed, three independent arguments):

1. **m_gamma and rho_M are not independent.** The canonical lattice relation is M_gamma^2 ∝ n_monopole
   (Anishetty et al / monopole-condensate literature). So where monopoles are dilute, the photon is light, and
   vice versa — they move together, not orthogonally. As a 2D separator they collapse toward one axis: they
   jointly separate {massless Coulomb} from {massive, monopole-screened} and cannot split the latter.

2. **m_gamma is massive in BOTH Confined and Higgs** (FS eqs. 3.6-3.8: Coulomb C(r)~1/r^4 power-law massless;
   Confined C(r)~exp(-4 ln(1/beta) r) massive; Higgs C(r)~exp(-m r) massive). So m_gamma = 0 is the unique
   Coulomb flag and cannot tell Confined from Higgs.

3. **rho_M is Higgs-screened to low in the deep-Higgs Z_q corner.** The charge-q condensate screens U(1)
   monopoles (measured: q=8, beta=0.8, rho_M 16.6@kappa=0 → 0.33@kappa=1.0). There the confinement is by Z_q
   electric flux, not monopoles, so rho_M stays low even though sigma_1 > 0. With m_gamma also massive there,
   {m_gamma, rho_M} read identically in deep-Higgs-Higgs and deep-Higgs-Z_q-confined → the confined phase
   artificially terminates at kappa~0.5. This is a genuine blind-order-parameter trap (same class as the
   SU(2)→2I case where L_link was blind and sigma_fund was the fix).

The deep-Higgs Z_q-confined region is REAL and paper-confirmed: Bowler fig.2/4 found the q=2 confinement line
persisting to kappa=inf → 0.4407, located by ENERGY HYSTERESIS; Damgaard-Heller traced the q=6 structure to
kappa=inf via the CJR Z_q transition lines. The screening selection rule (screened iff q|m, literally
src/u1/u1.hpp:135) is exact group theory, not a numerical artifact: for q≥2 charge-1 is never screened, so its
area law MUST persist into the deep-Higgs region. No monopole density (screened) + photon mass (massive in both)
can resolve a distinction carried by the unscreened charge-1 sector.

WHERE EACH OBSERVABLE WORKS / FAILS (no ambiguity):

- **m_gamma (structure factor):** WORKS as the unique massless-phase flag → isolates Coulomb (incl. the q≥5
  wedge). FAILS to distinguish Higgs from Confined (massive in both). HARD requirement: spatial extent L_s ≥ 16
  (at L=8 phat^2_min=0.586, at L=12 ≈0.27, both above the Higgs/screening scale → S(p) white → Coulomb and
  Higgs both read m^2≈0; at L=8 it is also wrong-signed in places). Use 16^3 × L_t geometric anisotropy. The
  temporal correlator is invalid — never plot it.

- **rho_M (monopole density):** WORKS as a confinement sub-classifier only in the strip kappa ≲ 0.5 (below the
  condensate) and beta ≳ 0.5 (above the DeGrand-Toussaint UV-noise floor). FAILS (screened low) for all kappa ≳ 0.5
  at every q → BLIND to deep-Higgs Z_q-confinement. DEMOTED from identifier to corroborator: its collapse while
  sigma_1 stays finite is the positive evidence that the residue is Z_q (flux-driven), not trivial.

- **sigma_1 (charge-1 Creutz string tension):** WORKS as the deep-Higgs Z_q-confined identifier — area law in
  BOTH Confined and deep-Higgs Z_q-confined (unscreened for q≥2), perimeter in Coulomb and ordinary Higgs. For
  q=1 it is 0 everywhere in the matter region (built-in screened control confirming the FS analytic connection).
  It is THE observable that does not terminate the confined phase at kappa~0.5. CAVEAT: "=0" is not literal —
  the finite-L perimeter+Coulombic residual chi(R,R) is ~0.01-0.06; use a quantitative gate (below), not an
  exact zero. Quotable sigma_1(beta) needs L ≥ 12 with the noise-guarded plateau; at L=8 only chi_1(2,2) is
  noise-reliable.

FINAL MINIMAL SET: IDENTIFY = {m_gamma, sigma_1, rho_M}, never any pair. sigma_1/|<P_1>| are REQUIRED for every
q≥2 (specifically on the deep-Higgs Z_q line and to certify the triple point). q=1 needs neither sigma_1 nor a
3-phase reading (2 phases, FS-connected). The Fredenhagen-Marcu order parameter rho_FM = G(R,T)/sqrt(W(R,2T))
(→0 in Coulomb, finite in both Higgs and Confined) is the gold-standard modern upgrade and the recommended
addition to make the Coulomb/Higgs cut at L ≤ 12 where m_gamma is dead — see §3.

---

## 3. LINE LOCATION (which observable locates each of the 3 lines + the triple point)

Reconciled with the 1980s recipe: Bowler-1981 and Damgaard-Heller-1989 drew EVERY line from (a) mean gauge &
link energies + their HYSTERESIS and specific-heat peaks, (b) extrapolation of the known pure-U(1) (beta_c≈1)
and pure-Z_q (CJR) endpoints inward, and (c) an external-field (Meissner/flux) probe to separate Coulomb from
Higgs. They did NOT use a photon-mass correlator. m_gamma and rho_M are modern refinements.

| Line | Primary locator | 1980s-paper backing | Where an ORDER PARAMETER is mandatory (susceptibility fails) |
|---|---|---|---|
| **Confinement/Coulomb beta-line** (small-kappa) | chi_plaq = Var(A)/V peak + cold/hot HYSTERESIS in <A>; anchor top at pure-U(1) beta_c≈1.01 | Bowler/D-H specific-heat peaks + hysteresis-loop width in <E> | — (gauge susceptibility scales; this line is robustly locatable at L≤12) |
| **Higgs kappa-line** (Coulomb/Higgs side) | m_gamma massless↔massive crossing at L_s ≥ 16; AT L ≤ 12: Fredenhagen-Marcu rho_FM →0 / external-field Meissner (penetrate=Coulomb / expel=Higgs) | D-H sect.3 + fig.5 external-field Meissner = the operational Coulomb/Higgs discriminator (done on 6^4) | **MANDATORY.** chi_link is a monotone q-blind CROSSOVER that does NOT scale (chi_link 4.80@L=8 = 4.73@L=12 — direct FSS). It LOCATES where matter energy turns over, never the phase boundary (Elitzur). Use rho_FM or the Meissner flux probe, not chi_link. |
| **Deep-Higgs Z_q line** (large-kappa, q≥2) | sigma_1 area↔perimeter crossing (+ |<P_1>|/chi_{P1} cross-check); anchor to kappa=inf pure-Z_q endpoint (q=2: 0.4407; CJR for q≥5) | Bowler q=2 line → 0.4407 via <E> hysteresis; D-H Z_q transition lines traced from CJR endpoints | **MANDATORY.** rho_M and m_gamma are blind here; only sigma_1/|<P_1>| (and rho_FM) see it. For q≥6 matter hysteresis is dead (gap=noise at L=8) — read the gap in the ORDER PARAMETER (sigma_1/|<P_1>|), not in <cos>_link, and use multicanonical-in-B to tunnel the first-order jump. |
| **q≥5 upper Z_q transition** (wedge→screened, beta~O(q^2)) | chi_plaq peak + cold/hot hysteresis in <A>, anchored to CJR pure-Z_q endpoint (beta_2~O(q^2)); m_gamma massless↔massive at L_s≥16 confirms | D-H fig.2/4 + CJR; located by gauge specific-heat/hysteresis, NOT a photon mass | gauge susceptibility suffices to LOCATE; m_gamma (L_s≥16) only confirms massless↔massive |
| **Triple point** (task 28) | intersection of (chi_plaq/sigma_1 drop) + (m_gamma or rho_FM Coulomb line) under FSS in L=8,12,16. q≤4: finite (beta_t, kappa_t). q≥5: recedes to kappa_t=inf (the wedge persists). | meeting of the gauge transition line and the Coulomb/Higgs line in both papers | Determine from the TWO measurable lines (chi_plaq+sigma_1, and rho_FM/Meissner); use m_gamma only as an L_s≥24 cross-check, not a defining line. |

Mandatory-order-parameter summary: the Higgs kappa-line and the deep-Higgs Z_q line CANNOT be located by any
susceptibility (Elitzur — no local Higgs order parameter; chi_link is a non-scaling crossover). They REQUIRE an
order parameter: sigma_1/|<P_1>| for the Z_q confinement cut, and m_gamma (L_s≥16) or Fredenhagen-Marcu /
external-field Meissner (L≤12) for the Coulomb cut. The two gauge lines (confinement/Coulomb beta-line and the
q≥5 upper Z_q transition) DO scale and are locatable by chi_plaq + hysteresis at L≤12.

---

## 4. IMPLEMENTATION (exact calls to wire into src/u1_frozen.cpp)

Frozen sampler is the correct vehicle: it has no modulus wall, matches Bowler/D-H fixed-length, and tunnels the
deep-Higgs gauge sector via the dB=0 Z_q heatbath (verified zq_flip ≈ 0.33 at kappa=3.0) — so its sigma_1 is
trustworthy where the HMC path (which cannot tunnel the freezing barrier) is not. `mode_point` currently emits
ONLY A, B, <plaq>, <cos>_link, chi_plaq, chi_link (confirmed at src/u1_frozen.cpp mode_point). The full IDENTIFY
suite is implemented and validated on the HMC scan path (src/u1_scan.cpp) and is a mechanical port — no new physics.

INCLUDES (order matters — monopole.hpp defines reduced_plaq_angle and sets the guard photon_structure.hpp respects):
```
#include "u1/monopole.hpp"           // FIRST. monopole_density<D>; defines reduced_plaq_angle<D>
#include "u1/photon_structure.hpp"   // photon_momenta, photon_structure_factor, photon_mass_fit (VALID m_gamma)
#include "measure/creutz_jack.hpp"   // CreutzJackknife, chi_diag_jack, plateau_sigma_jack, creutz_excl_reason, WGrid
#include "measure/observables.hpp"   // Stats (mean, binned_error, susceptibility(V))
```

PORT two free helpers VERBATIM from u1_scan.cpp into u1_frozen.cpp's anonymous namespace (self-contained on th/lat):
- `wilson_grids<D>(th, lat, q, Rmax, g1, gq)`  — u1_scan.cpp:52-90. Builds charge-1 grid g1[R][T]=<cos(loop)>
  AND charge-q screened-control grid gq[R][T]=<cos(q·loop)> in one OpenMP sweep.
- `polyakov_abs<D>(th, lat, n, tdir)`           — u1_scan.cpp:96-109. The REAL |<P_n>| modulus (answers the
  |P_1| question; do NOT use u1.hpp:155 `polyakov(...,m)` which returns <cos(line)>≈0).

IN `mode_point` measurement loop (after `s.sweep()`), on `s.th` / `s.lat` (kDim = 4):
```
rhoM.add( u1::monopole_density<kDim>(s.th, s.lat) );                    // rho_M
absP1.add( polyakov_abs<kDim>(s.th, s.lat, 1) );                        // |P_1|
wilson_grids<kDim>(s.th, s.lat, q, Rmax, w1blk, wqblk);                 // accumulate per-config grids
//   every n_block configs: jk1.add_block(w1blk); jkq.add_block(wqblk);
// L_s >= 16 build only: perConfig.push_back( u1::photon_structure_factor<kDim>(s.th, s.lat, mom) );
```
Build `mom = u1::photon_momenta<kDim>(s.lat)` ONCE before the loop. `rhoM`/`absP1` are `Stats`.

FINALIZE (copy u1_scan.cpp:287-307 verbatim — the reliableR noise gate + plateau):
```
auto reliableR = [&](const CreutzJackknife& jk, const WGrid& Wm, const WGrid& We, std::vector<int>& out){
  for (int R = 2; R <= Rmax; ++R) {
    if (!creutz_excl_reason(Wm, We, R, R, /*n_sigma=*/2.0).empty()) continue;
    const JackResult c = chi_diag_jack(jk, R);
    if (c.ok && c.value > 0.0 && c.value >= 2.0 * c.error) out.push_back(R);
  }
};
std::vector<int> R1, Rq; reliableR(jk1, jk1.mean_grid(), jk1.err_grid(), R1);
                          reliableR(jkq, jkq.mean_grid(), jkq.err_grid(), Rq);
const JackResult s1 = plateau_sigma_jack(jk1, R1);   // sigma_1 +/- err
const JackResult sq = plateau_sigma_jack(jkq, Rq);   // sigma_q +/- err
const Real absP1_mean = absP1.mean(), absP1_err = absP1.binned_error();
const Real chiP1 = absP1.susceptibility(Vspatial);   // Vspatial = lat.vol / L[tdir]
const Real rhoM_mean = rhoM.mean(), rhoM_err = rhoM.binned_error();
// L_s >= 16: PhotonMassFit fit = u1::photon_mass_fit<kDim>(perConfig, mom, nfit);  -> fit.m_gamma, fit.m2_err
```

NEW mode_point ARGS (optional trailing, like u1_scan.cpp): `Rmax` (default L/2, keep ≤ L/2 — loops wrap) and
`n_block` (default ~20). Emit one CSV row so the campaign analyzer parses it:
`beta,kappa,L,q, plaq±, <cos>_link±, A, B, chi_plaq, chi_link, rho_M±, sigma_1±, sigma_q±, |P_1|±, chi_P1 [, m_gamma± if L_s≥16]`.

NEW small observable: NONE strictly required — `polyakov_abs` already is the |<P_1>| order parameter. RECOMMENDED
addition for the Coulomb/Higgs cut at L ≤ 12 (where m_gamma is dead): the **Fredenhagen-Marcu** ratio
rho_FM(R,T) = G(R,T)/sqrt(W(R,2T)) (→0 in Coulomb, finite in Higgs AND Confined) — built from the SAME wilson_grids
plus the half-loop matter-line correlator G; ~40 lines, the gold-standard gauge-invariant Coulomb/Higgs separator.
Equivalent historical fallback: the D-H external-field Meissner probe (impose theta_ext = 2*pi*n_ex/L_perp^2,
measure flux penetration vs expulsion) — resolved on 6^4, so it works where m_gamma cannot.

GATES (turn the table into a classifier — no exact zeros):
- COULOMB:  m_gamma ≤ 2σ above 0 at L_s≥16  (OR rho_FM ≈ 0 at L≤12)  AND  sigma_1 < 2σ and sigma_1 < 0.1·sigma_1(confined-ref).
- CONFINED: sigma_1 ≥ 2σ above 0 (plateau) AND rho_M high (only checkable at kappa<0.5).
- DEEP-HIGGS Z_q-CONFINED (q≥2): sigma_1 ≥ 2σ above 0  AND  sigma_q < 2σ (screened control)  AND  rho_M low.
- HIGGS (Z_q-deconfined): sigma_1 < 2σ  AND  |<P_1>| > 0 (chi_{P1} below its peak)  AND  m_gamma massive (L_s≥16).
- q=1 control: sigma_1 < 2σ throughout the matter region → confirm FS analytic connection (no Z_q tongue, no triple point).

THE SCAN that renders all 3 phases unambiguous (per q ∈ {1,2,4,6,8}):
- WINDOW: beta from 0.2 to 2.5 (the Z_q-confined corner is at beta < 0.5; Z_2 endpoint 0.4407 — current grids
  stopping at beta=0.5 miss it), kappa from 0 to 3 (deep-Higgs needs kappa up to ~2-3; current grids stop at 1.2 miss it).
  For q≥5 also probe beta up to ~O(q^2) (the upper Z_q transition / wedge ceiling).
- VOLUMES: L = 8 (cheap survey, sigma_1 via chi_1(2,2) only, |<P_1>| sharp), L = 12 (quotable sigma_1 plateau),
  and a dedicated anisotropic 16^3 × L_t (L_t=6-8) run for m_gamma (the ONLY way to claim a Coulomb line).
- SAMPLER: frozen (u1_frozen mode_point) everywhere in the deep-Higgs corner. HOT+COLD starts and ≥2 seeds at
  every deep-kappa point; quote only where the brackets close. Use mode_muca (multicanonical-in-B) to tunnel the
  first-order Z_q jump for q≥6 where hysteresis in <cos>_link is dead.
- VERIFICATION GATE (run after wiring; the fingerprint must reproduce):
  (a) beta=2.5, kappa=0, q=2  (Coulomb)            → rho_M low, sigma_1≈perim, |<P_1>| large, m_gamma=0 @ L_s≥16.
  (b) beta=0.5, kappa=0, q=2  (confined)           → rho_M high, sigma_1>0, |<P_1>|≈0.
  (c) beta=0.9, kappa=2.0, q=4 (deep-Higgs Z_4)    → sigma_1>0 WITH sigma_q≈0 AND rho_M collapsed AND |<P_1>|≈0
       — this is the point the old {m_gamma,rho_M}-only map got WRONG (it terminated the confined phase at kappa~0.5).

---

## 5. HONEST RESIDUAL RISKS at L ≤ 12 — and the commitment

The protocol above is COMMITTED. These are the known limitations to manage, not reasons to relitigate:

1. **m_gamma is non-functional below L_s ~ 16.** At L=8/12 phat^2_min sits above the Higgs/screening scale, so
   the structure factor is white and both Coulomb and Higgs read m^2≈0 (at L=8 it is even wrong-signed in spots).
   MITIGATION: never claim a Coulomb line from L≤12 m_gamma; do the dedicated 16^3×L_t anisotropic run, and use
   Fredenhagen-Marcu / Meissner as the L≤12 Coulomb/Higgs separator. This is the single biggest gap.

2. **sigma_1 at L=8 is perimeter-saturated and can read falsely area-like** (short-distance self-energy the
   Creutz ratio cannot subtract). MITIGATION: quote sigma_1 only from the L≥12 noise-guarded plateau
   (reliableR + plateau_sigma_jack), always alongside sigma_q (screened control) and |<P_1>|, and use a
   quantitative gate (< 0.1·confined-reference), never an exact-zero test.

3. **HMC metastability poisons single-chain order parameters in the deep-Higgs q≥2 corner** (hot vs cold <plaq>
   differ by ~0.46 at the trapped points; sigma_1 non-monotone in kappa = sampling, not physics). MITIGATION:
   use the FROZEN sampler (tunnels via dB=0 Z_q heatbath), mandatory hot+cold + multi-seed bracketing, and
   multicanonical-in-B where brackets don't close. Treat the hot/cold gap itself as the first-order flag.

4. **rho_M is UV-noisy below beta~1**, exactly where the deep-Higgs Z_q tongue starts — so it cannot even serve
   as the low-kappa confined corroborator there. MITIGATION: rely on sigma_1/|<P_1>| in that corner; use rho_M
   only as a corroborator at beta≳0.5, kappa≲0.5.

5. **The deep-Higgs Z_q region and the q≥5 wedge have never been measured in the right window** (current grids
   stop at beta=0.5, kappa=1.2). MITIGATION: the extended window in §4 is non-optional; the central q≥2
   persistence claim must actually be sampled before any triple point or wedge is published.

Commitment: identify phases with {m_gamma (L_s≥16), sigma_1, rho_M} + |<P_1>| cross-check + sigma_q control;
locate the two gauge lines with chi_plaq + hysteresis, and the two matter cuts with sigma_1/|<P_1>| (Z_q line)
and Fredenhagen-Marcu / m_gamma@L_s≥16 (Coulomb line). sigma_1 is mandatory for every q≥2. This is the 1980s
recipe (energies + hysteresis + endpoint anchoring + an external Coulomb/Higgs probe) upgraded with the modern
order parameters, and it is settled.

---

### Memory actions (write-back)
- KEEP `u1-deep-higgs-zq-blindspot.md` — correct and load-bearing.
- EDIT `u1-coulomb-locator-photon-mass.md` — retract "m_gamma + rho_M separate all three phases"; replace with
  "m_gamma isolates Coulomb (unique massless, needs L_s≥16); rho_M is a confinement sub-classifier valid only at
  kappa≲0.5, beta≳0.5; the deep-Higgs Z_q tongue (q≥2) needs sigma_1/|<P_1>| (or Fredenhagen-Marcu)." Keep the
  rest (structure factor not temporal correlator; L_s≥16 wall; broken summary columns).
- `phase-diagram-plan-FS.md` 2-axis + m_gamma + sigma_1/|<P_1>| strategy = the correct synthesis; cross-reference this doc.
