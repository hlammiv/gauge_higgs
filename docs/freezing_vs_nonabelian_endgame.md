# Does freezing the radial mode cost us the nonabelian endgame?

**Question.** To beat the deep-Higgs first-order matter sampling wall in the U(1)+charge-q warm-up we are
considering freezing the radial mode (|phi|=1, lambda->inf). Does that freeze cost us anything we will *need*
for the real program: proving that a single large Higgs irrep R of SU(N) **selects** a specific discrete
nonabelian subgroup H (SU(2)->2T/2O/2I, SU(3)->Sigma(108)/...) so that kappa->inf IS the H gauge theory?

**Short answer.** No — *provided* "freeze" means what the code already does (per-site overall **scale**
|phi_x|=1, geodesic drift, **locking potential left ON**). That removes exactly the one gauge-invariant radial
direction and keeps the entire H-**selecting** shape orbit dynamical. The freeze that *would* destroy the
program — pinning the **shape** (the singlet-plane angle, the V_J ratios, the field onto a chosen singlet
vector or an H-orbit) — is a different, stronger constraint that no current code path performs and that we must
not build as a sampler. There are two real, bounded caveats: (i) a frozen sphere with **no** (or mis-tuned)
locking potential mis-selects to the principal abelian isotropy group, and (ii) radial-anchored Higgs
spectroscopy (the Observation-2 H-multiplet degeneracy fingerprint) is a lambda=inf casualty and must come from
a finite-lambda companion run.

---

## 1. CORE VERDICT — which directions are safe to freeze, which are load-bearing

Per site, a complex d-dim irrep field phi lives on C^d (2d real DOF). Decompose the gauge-invariant content:

- **The overall SCALE** |phi| — one real radial direction. This is the U(1)-analog modulus. **Safe to freeze.**
- **The SHAPE / directional sector** — the (2d-1)-sphere S^{2d-1}, on which the multi-invariant potential
  S_pot = sum_J f_J V_J (V_J = ||P_J(phi phi^dag)||^2, homogeneous degree-4, **direction-dependent**) selects
  **where on the orbit the VEV sits**, i.e. **which H**. **Load-bearing — must stay dynamical.**

H is the **little group** (stabilizer) of the VEV: G_phi = {g: D^sigma(g)phi = phi}. The VEV aligns along an
H-singlet of sigma|_H. H is a property of the *direction* of phi on the sphere, not its length. So:

| Operation | What it removes | H-selection | Verdict |
|---|---|---|---|
| **Overall-scale freeze** (|phi_x|=1, geodesic, potential ON) — *what GH_FROZEN does* | the radial ray only | untouched (shape orbit stays free) | **SAFE** |
| **Shape freeze** (pin theta / V_J ratios / a singlet vector / an H-orbit) | the H-selecting directions | inputs or mis-selects H | **FATAL** (mult>=2); question-begging always |
| **Scale freeze with NO / mis-tuned potential** | radial ray, but no selector left | relaxes to principal (abelian) isotropy, NOT H | **FATAL** |

### The decisive worked example: SU(2) spin-2 -> Q8 (multiplicity-2, "soft lock")

The Q8-singlet has **multiplicity 2** in spin-2, so the singlets span a 2D plane
phi(theta) = cos(theta) phi_1 + sin(theta) phi_2; every theta is Q8-invariant, but theta is **not** fixed by
group theory — it is selected dynamically by the quartic. Scanning theta (verified via the genuine 8-element
Q8 on spin-2): the gauge-orbit rank = number of broken SU(2) generators is **not constant**.
- generic theta -> **rank 3**: all three W's massive = the intended **discrete Q8** residual;
- special theta (theta = pi/2, 3pi/2, ...) -> **rank 2**: an **enhanced stabilizer** H' strictly containing Q8,
  with a **surviving continuous U(1) (a massless photon)**.

So the shape angle theta literally selects "discrete Q8" vs "Q8 + an extra unbroken photon," and which one wins
is fixed by the couplings f_J — exactly draft-1.tex:165. **A shape freeze pinning theta to the rank-2 point
yields a continuous residual and destroys the digitization claim.** This is the prompt's exact danger, made
quantitative. *But* an overall-scale freeze does **not** touch theta — it is an angle on the frozen sphere and
stays dynamical under geodesic drift — so |phi|=1 freezing is still safe even here, as long as theta is sampled.

### The canonical lockers: SU(2) spin-3 -> 2T (and SU(3) -> Sigma(108)) are multiplicity-1, RIGID

For the canonical lockers the H-singlet is **multiplicity 1**: SU(2) **2T at j=3**, 2O at j=4, 2I at j=6;
SU(3) **Sigma(108) at (2,2)**, Sigma(216) at (4,1), Sigma(648) at (3,3), Sigma(1080) at (6,0). A mult-1 singlet
is a **unique direction** (a Michel critical orbit of maximal stabilizer): dV/dphibar = c·phi for *every*
G-invariant potential, so H is fixed by **group theory alone** — the f_J only tune the radial gradient (mu^2)
and enforce stability/BFB. There is **no shape modulus to lose**. (Verified: 2T j=3 mult-1, orbit rank 3, full
discrete break, M_ab eigenvalues [4,4,4]; Sigma(108) (2,2) mult-1, M_ab = I.) For these, freezing the scale
carries **zero** selection risk — even a crude pin onto the singlet ray could not mis-select.

**Rule of thumb (the freeze-safety gate): mult-1 -> unconditionally shape-safe; mult>=2 -> still scale-freeze,
but keep theta sampled and verify the minimum lands on the max-orbit-rank discrete stratum.**

---

## 2. LEDGER — what we lose for the nonabelian goal

Columns: **SEL** = hurts H-selection? **ID** = hurts H-identification/proof? **MATCH** = hurts kappa->inf
matching to pure-H? Tags: NO / YES / partial.

### (a) Naive |phi|=1 freeze (scale-only **but with the wrong assumption that the sphere auto-selects H**)

| Item | SEL | ID | MATCH | Recoverable how |
|---|---|---|---|---|
| **Bare-quartic / no locking potential** on the frozen sphere | **YES** — relaxes to principal (abelian) isotropy, not H (research-program.md:16) | YES (wrong H identified) | YES (matches wrong theory) | **Keep S_pot = sum_J f_J V_J ON and copositive** on the sphere; this is the entire fix |
| Mis-tuned / non-copositive f_J (BFB is copositivity, not f_J>=0; "SU(2) BFB automatic" is REFUTED) | YES (mis-selects / unbounded) | YES | YES | Stay in the validated copositive cone (docs/locking_couplings.md) |
| **Global-minimum check** (Michel gives an extremum, not the global min; not yet verified) | partial — a wrong stratum could be the true global min | partial | partial | Must **sample** the sphere (never pin) to confirm the intended stratum wins |

### (b) Overall-scale-only freeze done right (GH_FROZEN, geodesic drift, potential ON, sphere sampled)

| Item | SEL | ID | MATCH | Recoverable how |
|---|---|---|---|---|
| H-selection (little group of the VEV) | **NO** | NO | NO | — (shape orbit stays dynamical) |
| Shape modulus theta in mult>=2 cases (Q8 j=2, 2T j=6, SU(3) (7,1)) | **NO** (theta stays free on the sphere) | NO | NO | Keep theta sampled; verify max-orbit-rank discrete stratum |
| Gauge-boson mass / unbroken-generator counting (n^a = Re[phi^dag T^a phi], links) | NO | **NO** — no radial derivative; needs only unit direction | NO | Runs frozen unchanged; cleaner (phi pre-normalized) |
| kappa->inf MATCHING to pure-H (sigma_fund, |P_fund|, rep-Polyakov, chi_plaq) | NO | NO | **NO — IMPROVED**: frozen removed the spurious radial-condensation kappa~0.3 ridge | This is the freeze's main payoff; the matching *wants* to be frozen |
| Elitzur-safe 2-axis classifier (chi_plaq/sigma_fund + L_link/chi_link) + screening table | NO | NO (all survive; L_link = clean directional hopping at \|phi\|=1) | NO | Runs frozen unchanged |
| **H-multiplet Higgs SPECTROSCOPY** (Obs-2 degeneracy fingerprint) | NO | **partial/YES** — radial mode gone, massive H-multiplet scalars stiffen (m_H->inf) | NO | **Finite-lambda companion run** (the draft already defers full spectroscopy); tangential/Goldstone correlators still survive frozen |
| Radial sigma mass; lambda-axis; **frozen-vs-unfrozen digitization-error** observable | NO | partial (secondary deliverables) | NO | Keep the unfrozen finite-lambda path alive; one control run per target H |

**Net:** the only genuine identification casualty of the correct scale-freeze is **radial Higgs spectroscopy**,
which we already planned to route to a finite-lambda companion. Everything load-bearing for **selection** and
**matching** survives, and matching is cleaner frozen.

---

## 3. THE RIGHT STRATEGY — what "freezing" should mean for the nonabelian scalar

**Definition of FREEZE (production rule):**
> Impose **one** per-site constraint |phi_x| = 1 (scale only). Keep the full multi-invariant locking potential
> sum_J f_J V_J **ON** and copositive. **Sample** the whole sphere S^{2d-1} — never pin the direction, the
> singlet-plane angle theta, the V_J ratios, a singlet vector, or an H-orbit.

This is exactly the existing `GH_FROZEN`/`frozen_phi` path: `project_pi_tangent()` removes only
Re(phi^dag pi)phi (the radial momentum), `geodesic_drift()` is the exact reversible great-circle flow on the
sphere, and `kick()` evaluates the **full** scalar force including the MultiInvariantPotential gradient *before*
projecting — so the potential acts on the shape sector while only the radial component is dropped
(src/hmc/gauge_higgs_hmc.hpp:63-103,133,144; src/action/scalar_invariants.hpp:60-73). Because V_J are
homogeneous degree-4, the minimizing **direction** and the entire validated f_J copositive cone transfer to the
frozen theory unchanged — **no re-derivation of locking couplings is needed**.

**An explicit H-orbit constraint is NOT a sampler.** Pinning phi onto the H-singlet plane or the orbit G·v_H
inputs H by construction, begs the very question "R + potential SELECTS H," hides a wrong-stratum global
minimum, and erases the soft-lock modulus that distinguishes Q8 from a rigid lock. It is legitimate **only** as
a kappa->inf *analysis/control* limit (links freezing onto the stabilizer = residual H, the matching target),
**never** as the sampler used to *prove* selection.

**Implication for task #13 (the "RATTLE" item).** The closed-form geodesic |phi|=1 integrator is the **right**
and already-shipped nonabelian scalar update — it is exact, reversible, |phi|^2=1 to machine precision,
<e^-dH> -> 1, and is FD-verified on the production GeneralRep spin-3 + MultiInvariantPotential path at kappa=8
and 64. Two things to internalize:
- **Single-|phi|=1 RATTLE is the U(1) special case, NOT a nonabelian-specific tool.** In U(1) the field is one
  complex number, so |phi|=1 removes the *only* modulus. In a large irrep the same constraint removes one of
  2d directions and leaves the (2d-1)-shape sector — that is *why* it stays safe, not because the modulus is
  unique. State the bridge precisely: harmless because **only the scale is frozen, the H-selecting shape stays
  alive**, not because "the modulus is unique."
- **Do NOT build a multi-constraint / invariant-ratio RATTLE.** Fixing several V_J values would freeze the
  shape space (the load-bearing directions), inherit a curved phi-dependent constraint manifold with iterative
  solver cost, and mis-select H for every mult>=2 rep. It is redundant at best, fatal at worst.

**Sampling, by sector (freezing does not cure the real wall).** Freezing helps the **matter** sector (it
removed the radial-condensation artifact) but is **neutral** on the actual deep-kappa wall, which is the
**gauge-sector first-order freezing** of the plaquette (hot/cold L_link agrees while the plaquette stays
bistable). Stack tools by sector:
- *matter*: frozen geodesic-sphere HMC (default);
- *gauge*: compose **GH_GUPD center-flip** sweeps between trajectories — exact (dS_H == 0) for every locking
  rep because all are **center-blind** (integer-j SU(2); p+2q == 0 mod 3 SU(3)); KP heatbath leg is N==2 only,
  so SU(3) needs Cabibbo-Marinari, with the N-general Metropolis leg as interim fallback;
- *deepest kappa*: **multicanonical / LLR** in the plaquette (machinery exists, build/u1_llr).

**Compute note.** Freezing is cost-neutral on the hot path; the dominant MD cost is the locking-potential
gradient (cured ~19x by the Asuper superoperator hoist) and the rep apply (cured by fast_D = exp(i w·T_R)).
Lead SU(3) with Sigma(108) (2,2), d=27 (Asuper 729x729, cheap); flag a d=64 Asuper feasibility check before
committing to Sigma(648).

---

## 4. BOTTOM LINE — is the U(1) frozen warm-up still the right first step?

**Yes.** Freeze the U(1)+charge-q warm-up: it is the genuinely harmless single-modulus case (one complex field
per site; H = Z_q is fixed by the **charge / center (N-ality)** in the e^{iq theta} hopping, independent of
|phi|), it beats the matter sampling wall, and the same geodesic-sphere machinery transfers directly to the
nonabelian scalar. Then freeze the **SU(2)->2T (j=3) headline** and the **kappa->inf matching** as well — all
mult-1, all shape-safe.

**What NOT to over-learn from U(1):**
1. **"One modulus fixed by charge" does NOT transfer.** The nonabelian center/N-ality only separates
   center-blind from center-active reps; it does **not** pick the specific H among 2T/2O/2I (all center-blind).
   H is selected by **where the VEV sits in the irrep**, i.e. by the shape sector + the f_J — not by a charge.
2. **"Freezing removes nothing" is false in general.** It removes nothing **for selection** only because the
   shape stays dynamical **and** the locking potential stays on. A bare-quartic frozen sphere mis-selects to
   the principal abelian group. Freezing is harmless **conditional** on keeping the H-selecting potential alive,
   tuned (copositive), and the sphere **sampled**.
3. **Shape moduli are load-bearing whenever multiplicity >= 2.** Treat Q8 spin-2 as the deliberate canary
   (never the headline locker); never pin theta or the V_J.
4. **Keep one unfrozen finite-lambda control per target H** — for Higgs spectroscopy (Obs-2 degeneracy), the
   lambda-axis, and the frozen-vs-unfrozen digitization-error observable.

**Operational one-liner:** *Freeze |phi_x|=1 only, never the direction; keep sum_J f_J V_J on and copositive;
gate every irrep with a singlet-multiplicity check; sample the sphere; route radial spectroscopy to a
finite-lambda companion.*

---

*Anchors:* draft-1.tex:111,127,165,270,274-276,580-588,767,915-975; research-program.md:16; locking_couplings.md:17-18,25,47-58;
phase-diagram-plan-FS.md:14-16; su2-discrete-controls.md; su2-discrete-three-phases.md; discrete-betac-hmc-cannot-tunnel.md;
gauge-boson-mass-fms-gauge-variance.md; src/hmc/gauge_higgs_hmc.hpp:63-103,133,144; src/action/scalar_invariants.hpp:60-73;
src/measure/gaugeboson_op.hpp; src/screening.cpp; src/hmc/gauge_link_updater.hpp; plus the Q8 spin-2 theta-scan
(rank-3 Q8 vs rank-2 enhanced-stabilizer+photon) and 2T j=3 mult-1/rank-3 verification.
