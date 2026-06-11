# How q=2 / q=6 were extracted 40 years ago — and what it means for our wall

**Decision brief for the lead.** Two papers mapped the U(1)+charge-q Higgs phase diagram (q=1,2 frozen; q=2,6 variable) on 6⁴–8⁴ lattices with 1980s hardware, while our HMC + 2D-LLR + multicanonical stack stalls on the deep-Higgs matter axis. This explains exactly why, separates "located the lines" from "determined the order," and gives a verdict on freezing |φ|=1.

Anchor papers:
- **Bowler, Pawley, Pendleton, Wallace, Thomas**, *Phase Diagrams of U(1) Lattice Higgs Models*, Phys. Lett. B 104B (1981) 481 — `literature/1-s2.0-0370269381905190-main.pdf`. **Fixed-length**, q=1,2, 8⁴, DAP Metropolis.
- **Damgaard & Heller**, *Vortices and the Phase Structure of the Multiply-Charged U(1) Higgs Model*, Nucl. Phys. B324 (1989) 532 — `literature/1-s2.0-0550321389904793-main.pdf`. **Variable-length**, Q=2,6, λ=1.0 and 0.01, 6⁴/8⁴, Metropolis.
- **Fradkin & Shenker**, Phys. Rev. D19 (1979) 3682 — actual file `literature/phase_diag_LF.pdf` (the `PhysRevD.19.3715.pdf` in the folder is Horn-Weinstein-Yankielowicz, not FS). **Frozen-length** analyticity theorem.

Our action (`src/u1/u1.hpp:65`): `S_H = Σ_x[|φ|² + λ(|φ|²−1)²] − κ Σ 2Re[conj(φ)e^{iqθ}φ′]`, dynamical modulus, λ=0.5 default (`src/u1_scan.cpp:350`).

---

## Q1 — How did they succeed where our HMC + LLR + multicanonical fails?

**Short answer: they never tried to do the thing that is breaking us.** They (a) used a sampler class (local Metropolis with hand-built discrete moves) that can teleport across a first-order barrier, (b) located lines by hysteresis / thermal cycling on tiny lattices, NOT by tunnelling a coexistence region, and (c) fixed the *order* of the transition from exact analytic limits, not from finite-size scaling. We are failing partly on physics (the dynamical-modulus barrier is genuinely worse) and **partly because we hold a strictly higher standard** (FSS-grade order determination by directly sampling the coexistence) than either paper attempted.

### The concrete recipe

| Ingredient | Bowler 1981 | Damgaard-Heller 1989 | Us |
|---|---|---|---|
| Modulus | **Frozen** \|φ\|=1 (one angle χ/site) | Variable length ρ·e^{iχ} | Variable length (dynamical ρ) |
| Update | Local single-angle **Metropolis** | Local **Metropolis** (modulus+phase) | **HMC** (MD trajectories) |
| Gauge trick | — | **Z_Q link jump** U→e^{±2πi/Q}U | none equivalent |
| Topology trick | — | **Global one-flux-unit** Metropolis move | none (HMC conserves sector) |
| Volume | 8⁴ only | 6⁴ (primary), 8⁴ | 4⁴–12⁴, pushing larger |
| Statistics | Δβ≈0.025, 32–100 sweeps/pt | Δκ<0.05, 50–100 iters/pt | long; FSS-grade target |
| Line location | **Hysteresis loop width** in ⟨S⟩/⟨E⟩ | **Thermal cycling**, jump in ⟨ρ²⟩ | trying to **tunnel** coexistence |
| Order | hysteresis persistence (admittedly crude) | **exact analytic limits** | trying FSS/LLR/muca |
| FSS? | **No** (disclaimed) | **No** (disclaimed) | yes (the hard part) |

### Why their samplers crossed the barrier and ours can't
1. **Metropolis can be handed a barrier-crossing proposal; HMC cannot.** A single tailored Metropolis proposal can land *exactly* on the far side of a first-order well with O(1) acceptance. HMC integrates a continuous symplectic trajectory — it is barrier-bound and sector-conserving by construction. This is precisely our established wall (local HMC can't tunnel latent-heat×volume).
2. **The Z_Q link-jump move (D-H, p.535 footnote).** U→e^{±2πi/Q}U lands a link directly on the neighbouring Z_Q minimum that the κ→∞ action selects. Because `e^{iqθ}` is *exactly invariant* under θ→θ+2π/q, the extensive κB term **cannot veto the move** — only the plaquette term can. They call it "advantageous, especially in the large κ region." We have *no analogue in the U(1) HMC* (our SU(2) GH_GUPD is the same idea; D-H validate it for U(1)+q).
3. **Global one-flux-unit Metropolis (D-H, p.539).** Accept/reject of an entire configuration changing total flux by 2π — a non-local collective move between metastable topological sectors. HMC's continuous flow fundamentally cannot make this jump.

### They LOCATED the lines ≠ they DETERMINED the order (be honest about both)

**Location (both papers): hysteresis / thermal cycling — NOT tunnelling.**
- Bowler: sweep β up then down (Δβ≈0.025, 32–100 sweeps/step), read transition off the **width of the hysteresis loop** in ⟨S⟩/⟨E⟩; "thick lines" on Figs. 1–2 are loop widths. They ran cold-start and hot-start sweeps and *never crossed the barrier in equilibrium*.
- D-H: thermal cycling in κ (Δκ<0.05, 50–100 iters/pt), read jump in ⟨ρ²⟩=⟨|φ|²⟩ (p.536). Located the β=0 line vs an exact V_eff (eq.4) — "agreement almost perfect."
- Both **explicitly disclaim finite-size scaling**. D-H: "we have not checked in detail the dependence on the lattice size." Bowler: single 8⁴, "can say little about the detailed nature of the phase transitions."

**Order (both papers): mostly ANALYTIC, not from MC.**
- D-H pinned order from exact corners: β=0 gauge-integrated V_eff(ρ) double-well → **first order for λ≤0.13** (eq.4); β=∞ → O(2) φ⁴ → **second order**; κ=∞ → imported pure Z_Q transitions (first order Q≤4; two transitions Q≥5, Creutz-Jacobs-Rebbi). MC only *connected and located* lines.
- Bowler inferred order from hysteresis persistence: strong first order at q=2, β_P=0.5, β_L=4.0 (Fig. 4, stably separated heating/cooling branches); "only weakly [first order]" near the low-β_L (triple-point) end. Honest caveat: "significantly better statistics would be required" to find where order changes.

**Net for us:** the barrier they hit when they *did* try to tunnel was real and beat them too — D-H's Q=6 deep-Higgs flux-sector metastability "exceeded our simulation time" (p.544 footnote). They simply **routed around it**: location by cycling + jump-spotting, order by analytic limits. We are demanding the one measurement (sample both coexisting phases → FSS the order) that *neither paper performed*. Part of our "failure" is a self-imposed higher bar.

---

## Q2 — How important is the fixed radial coupling?

**Verdict: it is the single highest-leverage change available. Freezing |φ|=1 dissolves the half of our wall that actually kills us (the radial slow mode + the unbounded, radially-amplified B), converts the remaining genuine first-order barrier into the textbook bounded-domain multicanonical problem the field has solved since the 1980s, and is fully physically legitimate for the triple-point / κ→∞ goal — arguably *more* faithful, because Fradkin-Shenker is a frozen-length theorem.**

### Why the radial mode is the mechanism of our wall
Our potential is `|φ|² + λ(|φ|²−1)²`, λ=0.5. Mean-field radial minimum with ordered hopping: ρ² = 1 + (2κDc−1)/(2λ) → ρ²≈1.9 at κ=0.3, 3.8 at κ=0.5, ~8 at κ=1 (D=4). Consequences, each tied to a documented failure:
- **B is unbounded and radially amplified.** B/link ~ 2ρ²c; angular alignment c↑ deepens ρ²↑ which multiplies κρ²c↑ — positive feedback. Latent ΔB is ~2–3× the frozen value, with an open-ended high-B tail (measured high-B branch 2335→6000 vs frozen ceiling 2·vol·D=2048 at L=4). Every method chased this tail: LLR windows to B~6000, WL bins over an unbounded range.
- **The radial condensate is the documented slow mode.** The LLR root cause ("the φ condensate does not reach equilibrium per RM iteration"; "|φ| modulus not equilibrated within the box"; 43/50 high-B cells failed at acceptance <0.3) is a *radial* statement: |φ| is a soft, extensive collective coordinate dragged by κ while HMC's step is capped by stiff angular/gauge forces — the worst case for MD.
- **Muca-in-force blew up on exactly this coordinate.** Bias g′(B) acting as κ_eff inside the integrator on an unbounded runaway B → κ_eff~−134 freeze (this is the `keff()` in-force path at `src/u1/u1.hpp:201`).

### What freezing removes vs keeps (honest)
- **Removes (the half that kills us):** the radial slow mode (→ LLR acceptance collapse gone), the unbounded radially-amplified B (→ muca/WL instability gone, no tail to chase), and HMC's step-size problem (matter → one compact angle/site with an *exactly samplable* conditional p(χ)∝exp(a·cos(χ−α)), a=2κ|Σ_neighbors|, exact at **any** κ → von Mises heatbath + rejection-free overrelaxation χ→2α−χ).
- **Keeps (genuine, irreducible):** the first-order character. κ=∞ terminates on the pure Z_q gauge transition (strongly first order for q≤4; Bowler Fig. 4 shows it strongly first order *frozen*). No local update — heatbath, overrelaxation, cluster — beats latent-heat×volume coexistence. **You still need multicanonical-in-B or LLR-in-B to cross it.** BUT now on a **bounded** domain |B|≤2·vol·D with finite bins and a per-move ΔB that is O(1)-bounded — the exact Berg-Neuhaus / Langfeld-Lucini configuration these methods were *invented and proven* in, and exactly what our own failure ledger already prescribed ("matter Metropolis update, muca bias in the accept, no integrator to destabilize").

### Physical legitimacy (Fradkin-Shenker topology argument)
Fully legitimate, and the frozen model **is** the FS model:
- FS abstract + Sec. I.B (`phase_diag_LF.pdf`): "we freeze out the radial [mode] … *without, we feel, throwing away any important physics*"; the Osterwalder-Seiler connectivity proof is "specialized to the case of fixed-length Higgs." Checking the triple point / connectivity against the frozen model is checking the theorem against its own assumptions.
- **κ→∞ digitization endpoint is λ-independent.** D-H's unitary-gauge argument (p.534): links freeze onto Z_Q for *any* λ; λ enters only through a renormalized ρ̄² (their eq.2). The residual discrete-group physics we want survives identically.
- **D-H demonstrated numerically** that λ=1.0 variable-length reproduces the frozen topology for *both* Q=2 and Q=6 ("qualitatively similar to the known λ→∞ case"). Our λ=0.5 is far above their λ≲0.13 threshold where a genuinely new radial line appears.
- **What freezing discards** (none in scope for #28/#26): the small-λ radial first-order line, the radial/Higgs-mass spectrum, and quantitative λ-dependence of line positions. The κ-axis is reparametrized (κ_frozen ≈ κ⟨ρ²⟩), so triple-point *coordinates* shift, but its **existence, order, and q-dependence are invariant**. Frozen is also the direct apples-to-apples comparison with both anchor papers.
- **FS topology we are testing:** q=1 → Higgs and confinement *analytically connected* (no boundary for β_P≲0.5, β_L≳1.0 — Bowler confirmed). q≥2 → fundamental-Wilson-loop confinement line *persists to κ=∞* (β=∞ reduces to a nontrivial Z_q gauge theory), giving the three-phase confined/Coulomb/Higgs structure and the triple point. The deep-Higgs confined phase is carried by a **gauge-sector** order parameter (charge-1 Wilson loop σ₁ / Polyakov |P₁|) that does **not** require resolving the matter B-jump — consistent with our `u1-deep-higgs-zq-blindspot` memory.

---

## What we should change — top moves (highest leverage first)

1. **Build a frozen-PHASE matter sector (task #13), NOT RATTLE-HMC.** *Effort: medium (~1–2 days).* Add site angles χ_x, B=2Σcos(qθ+χ′−χ); matter update = exact **von Mises heatbath** (Best-Fisher, a=2κ|resultant of gauge-rotated neighbours|, exact at any κ) + 3–5 **rejection-free overrelaxation** sweeps per heatbath sweep. *Reason:* deletes the radial slow mode, bounds B, removes HMC's step-size wall — the LLR acceptance collapse and the WL/muca blowup are gone *by construction*. **Do NOT use RATTLE-constrained HMC** — it keeps the force-based integrator and inherits exactly the failure mode we are escaping (and the U(1) frozen scalar is literally one angle/site, so the HMC machinery buys nothing).

2. **Add the Z_q gauge link-jump move θ→θ±2π/q (D-H, p.535).** *Effort: low (~hours; wire into BOTH frozen and the existing dynamical code).* *Reason:* the hopping term is *exactly invariant*, so the extensive κB term cannot veto it; it is the abelian twin of our verified SU(2) GH_GUPD center flip and decorrelates the deep-Higgs Z_q sectors (directly addresses the q≥4 residual-Z_q branch trapping in the L=8 tempered runs). Action-exact in the dynamical model too — low-risk free win regardless of #13.

3. **Adopt hysteresis / thermal-cycling for LINE MAPPING; reserve LLR/FSS for the ORDER question only.** *Effort: low (scripting on top of #1).* *Reason:* with frozen spin + heatbath+OR both branches equilibrate cleanly at any κ, so Bowler-style hot/cold bracketing is a *legitimate cheap locator* again — and the deep-κ line position is independently pinned by the analytic Z_q map (Bowler eq.5: β_eff=β_P(1−1/2β_L)→Z₂ self-dual 0.4407; D-H eq.2 generalizes to any Q). This matches the 1980s standard that succeeded; stop spending tunnelling effort just to *find* lines.

4. **Move the multicanonical weight g(B) into the ACCEPT step of the local frozen moves (bias-in-accept), drop the in-force path.** *Effort: low–medium (reuse `src/u1/u1_mucab.hpp`; bypass `keff()` in-force at `u1.hpp:201`).* Optionally re-point `build/u1_llr`'s matter sampler at the heatbath/Metropolis χ-update (per-move ΔB bounded → high-B acceptance collapse gone). *Reason:* this is the Berg-Neuhaus-proven configuration on a now-bounded domain — the genuine first-order barrier (which freezing does NOT remove) becomes the standard, routinely-crossed muca-in-B problem at L=4–12. Start with muca-in-accept; our own ledger already converged on it. **Skip embedded-Wolff clusters** — gauge frustration degrades them and clusters don't beat first-order coexistence.

**Validation ladder (all from the anchor papers):** pure-gauge β_c≈1.01 at κ=0; Bowler q=1 = no Higgs/confinement line for β_P≲0.5, β_L≳1.0 (FS connectivity control); Bowler q=2 deep-κ line → Z₂ self-dual β_P=0.4407 via the analytic map; D-H Q=6 intermediate massless phase reaching the κ=∞ edge (the q≥5 Coulomb-wedge headline); finally **one dynamical-λ=0.5 control slice** mapped by κ⟨ρ²⟩ to demonstrate frozen↔dynamical topology agreement for the paper.

---

## ⚠️ Does this REROUTE the separately-running large-κ sampling plan? — YES.

The current plan attacks the deep-Higgs matter axis in the **dynamical-modulus** model with **HMC-based** 2D-LLR / multicanonical-in-force / tempering. Every one of those failed on the *radial* coordinate (slow mode, unbounded B, κ_eff blowup, ladder variance). The recommendation **redirects the large-κ campaign to the frozen-phase model with local heatbath+overrelaxation**, demoting HMC for the matter sector entirely on that axis. Concretely:
- **Stop** investing in dynamical-modulus LLR matter-windowing and muca-in-force for the deep-Higgs line — they are fighting a coordinate that freezing deletes.
- **Keep** the dynamical-modulus runs only for: the validated gauge/confinement axis (β_c≈0.98–1.0), and **one** λ=0.5 topology-control slice for the paper.
- The frozen campaign is **not** a detour — it is the apples-to-apples match to Bowler/D-H/FS and the configuration in which the residual genuine barrier is a solved problem. The triple-point coordinates will shift (κ_frozen≈κ⟨ρ²⟩); existence, order, and q-dependence will not.
