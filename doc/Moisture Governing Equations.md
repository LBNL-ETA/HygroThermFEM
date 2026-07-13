# Moisture governing equations — reconstructed from the assembly code

> **Status: validated against the original THERMM documents (2026-07-10).** This document
> reconstructs the continuous PDE, its weak form, and the discrete system **directly from the
> element-assembly code** so the equations can be audited. The reconstruction was cross-checked
> against the reference documents *THERMM: Theoretical Model* (rev. 03/20/2020) and *THERMM:
> Numerical Implementation* (rev. 5/28/2019): the PDE (§2), weak form (§3), coefficient
> formulas (p_sat, R_v, D_phi = xi*D_l) and the moisture BC are all confirmed. D1 turned out
> to be an erratum in the Numerical Implementation document itself (its appendix §2.14.3 drops
> `phi div(delta grad c_sat)` between eq. 515 and the discretization); the thermal-side DpDu
> terms are advective per Theoretical eq. 424 and are **correct as written** (earlier suspicion
> retracted). Full cross-check, errata and decision log:
> `D:\Documents\HygroThermFEM Solver\moisture_solver_cross_check.typ`.

Reconstructed from: `src/elements/Element2D.cxx`, `src/elements/Element2D.hxx`,
`src/elements/Elements2D.cxx`, `src/math/Functions.cxx`, `src/materials/Material.cxx`,
`src/domain/Domain.cxx`.

---

## 1. Primary variable and model class

The moisture domain solves for **relative humidity** `φ` (`Variable::humidity`), stored per
node and clamped to `[0, 1]` after every solve (`MoistureDomain::postProcess`,
`MoistureDomain.cxx:230`). This is a **single-variable, RH-driven moisture-balance model** of
the Künzel/WUFI family: moisture storage is given by a sorption curve `w(φ)` and transport is
the sum of vapor diffusion and liquid (capillary) transport.

Symbols and code sources:

| Symbol | Meaning | Units | Source |
|---|---|---|---|
| `φ` | relative humidity (primary variable) | – | `Variable::humidity` |
| `w(φ)` | total moisture (sorption) storage function | kg/m³ | `sorptionCurve`, `Material.cxx:402` |
| `ξ = dw/dφ` | moisture storage capacity | kg/m³ | `TabularDerivativeSmooth`, `Element2D.cxx:590` |
| `μ` | vapor diffusion resistance factor | – | `diffusionResistanceFactor` |
| `δ = D_v/μ` | effective vapor diffusion coeff (concentration form) | m²/s | `Element2D.cxx:560`, `D_v = 2.5e-5` |
| `p_sat(T)` | saturation vapor pressure | Pa | `vaporPressureAtTemperature`, `Functions.cxx:49` |
| `c_sat(T)` | saturation vapor **concentration** = `p_sat/(R_v T)` | kg/m³ | `saturationConcentrationAtTemperature`, `Functions.cxx:56` |
| `D_l(w)` | liquid transport coefficient | m²/s | `LiquidTransportationCurve`, `Functions.cxx:378` |
| `R_v` | water-vapour gas constant (461.4) | J/(kg·K) | `Functions.cxx:59` |

Note the model works with saturation **concentration** `c_sat = p_sat/(R_v T)` (kg/m³), not
saturation pressure, so the vapor flux comes out directly in kg/(m²·s). This is dimensionally
equivalent to the pressure-based Künzel form via `δ_p ∇(φ p_sat) = (δ/R_v T)∇(φ p_sat) ≈
δ ∇(φ c_sat)` as long as the `∇(1/T)` contribution is treated consistently — which is exactly
where [Defect #1](#d1-inconsistent-weak-form-of-the-vapor-term-equation-level) lives.

---

## 2. Continuous PDE

Moisture mass balance (per unit volume):

```
ξ(φ) ∂φ/∂t  =  −∇·( g_v + g_l )
```

with the two mass fluxes exactly as the code builds them:

```
vapor:   g_v = −δ ∇(φ · c_sat(T))              (Element2D.cxx:558-571)
liquid:  g_l = −D_l(w) · (dw/dφ) ∇φ = −D_l ∇w  (Element2D.cxx:576-585)
```

so

```
ξ(φ) ∂φ/∂t  =  ∇·( δ ∇(φ c_sat) )  +  ∇·( D_l (dw/dφ) ∇φ )
```

The vapor driving force `∇(φ c_sat)` expands into a moisture-gradient part and a
**temperature-gradient part** (because `c_sat = c_sat(T)`):

```
∇(φ c_sat) = c_sat ∇φ  +  φ ∇c_sat
                          └── nonzero only where T varies (coupled cases)
```

The liquid coefficient `D_l` is evaluated from the material's own water content
`w = w(φ)` (`LiquidTransportationCurve::value`, `Functions.cxx:390`), not the node-averaged
value, to keep `K·φ_uniform = 0`.

---

## 3. Weak (Galerkin) form

Test function `ψ_i`, integrate over element `Ω`, integrate the flux divergence by parts once:

```
∫ ψ_i ξ ∂φ/∂t dΩ
   = −∫ ∇ψ_i · δ ∇(φ c_sat) dΩ  −∫ ∇ψ_i · D_l(dw/dφ) ∇φ dΩ
     + ∮ ψ_i (g_v + g_l)·n dΓ        (boundary → moisture BCs)
```

Expanding the vapor term with `∇(φ c_sat) = c_sat∇φ + φ∇c_sat`, the **consistent** interior
contribution is:

```
−∫ δ c_sat ∇ψ_i·∇φ dΩ        (A)  moisture-gradient vapor diffusion
−∫ δ φ    ∇ψ_i·∇c_sat dΩ     (B)  temperature-gradient vapor coupling
−∫ D_l(dw/dφ) ∇ψ_i·∇φ dΩ     (C)  liquid transport
```

Every term is integrated by parts once, so **every term differentiates the test function**
`∇ψ_i`. That is the consistency requirement a single-integration-by-parts Galerkin scheme
must satisfy. Hold onto that for §6.

---

## 4. Discrete system (backward Euler)

Per-element matrices (`IElementLinear2D`, integrator classes in `Element2D.cxx`):

| Code call | Integrator | Assembles | Continuous term |
|---|---|---|---|
| `DDu(k)` | `QLEDDuIntegrator2D` | `K_ij = ∫ k ∇ψ_i·∇ψ_j` | ✓ weak diffusion (test differentiated) |
| `DpDu(k, p)` | `QLEDpDuIntegrator2D` | `C_ij = ∫ k ψ_i (∇ψ_j·∇p)` | see §6 — **test NOT differentiated** |
| `Cap(c)` | `QLECapacitanceIntegrator2D` | `M_ij = ∫ c ψ_i ψ_j` | capacity `ξ ∂φ/∂t` |

Coefficients are applied as the edge-average `½(k_i + k_j)` at assembly
(`IQLEIntegrator2D::integrate`, `Element2D.cxx:23-44`).

Global assembly and time stepping (`Elements2D.cxx`, `Domain.cxx:225-249`):

```
LHS = K_DDu + C_DpDu + M/Δt              (conductanceMatrix sums DDu+DpDu — Elements2D.cxx:40;
                                          lumped mass on diagonal — getLumpedMass:115)
RHS = (M/Δt)·φⁿ + R_boundary + R_element
solve  LHS · φⁿ⁺¹ = RHS                   (backward Euler, 1st order, lumped/very diffusive)
```

The moisture element populates these as:

```
DDu( δ · c_sat )                 → term (A)   Element2D.cxx:565
DpDu( δ , c_sat )                → term (B)*  Element2D.cxx:567   (*see Defect #1)
DDu( D_l · dw/dφ )               → term (C)   Element2D.cxx:584
Cap( dw/dφ )                     → capacity   Element2D.cxx:594
```

The nonlinearity (all coefficients depend on `φ` and `T`) is resolved by the modified
Newton–Raphson loop in `Domain.cxx` (see `doc/IterationMethods.md`).

---

## 5. Coupling and the thermal side

`c_sat(T)` and `D_l(w)` make moisture depend on the thermal solution; the thermal element in
turn carries latent-heat / vapor-enthalpy terms that depend on `φ` (`Element2D.cxx:427-540`).
The two domains are advanced by a **single-pass staggered scheme** with adaptive sub-stepping
(`MultiDomain::transient`, `MultiDomain.cxx:45-207`) — fields are exchanged once per sub-step,
with no coupling (Picard) iteration to convergence within a step.

---

## 6. Defects found

### D1. Inconsistent weak form of the vapor term (equation-level)
The two halves of the single divergence `∇·(δ∇(φ c_sat))` are discretized with **different
integration-by-parts treatments**:

- Term (A) `δc_sat∇φ` → `DDu` → `∫ δc_sat ∇ψ_i·∇ψ_j` — integrated by parts, **test
  differentiated**. Correct.
- Term (B) `δφ∇c_sat` → `DpDu(δ, c_sat)` → `∫ δ ψ_i (∇ψ_j·∇c_sat)` — **test not
  differentiated** (`Element2D.cxx:105-118`). This is the Galerkin form of the *strong* term
  `δ ∇φ·∇c_sat`, i.e. it keeps only the first half of `∇·(δφ∇c_sat) = δ∇φ·∇c_sat +
  δφ∇²c_sat` and **silently drops `δφ∇²c_sat`**.

Consequences:
1. Mixing a by-parts term (A) with a not-by-parts term (B) for one divergence is
   inconsistent; the assembled operator is **not** the Galerkin discretization of
   `∇·(δ∇(φc_sat))`.
2. The consistent form of (B) is `−∫ δφ ∇ψ_i·∇c_sat` (test differentiated, transpose of what
   the code builds) — structurally different from `∫ δ ψ_i(∇ψ_j·∇c_sat)`.
3. The whole term is **temperature-gradient driven**: it is identically zero under uniform
   `T` and grows with `∇T`. This matches the reported symptom — *change the temperature
   (or a coupled parameter) and the moisture result swings*.

**Fix direction:** assemble the vapor flux as one consistent weak term
`−∫ δ ∇ψ_i·∇(φ c_sat)`, i.e. either put term (B) through the by-parts path (differentiate
`ψ_i`, keep `φ` from the previous iterate as a coefficient), or move it to a consistent RHS
source. Quantify the current error against the Python 1D reference with an imposed `∇T`
([Phase B, scenario 2](#)).

> **Document cross-check (2026-07-10):** the Theoretical Model (eq. 36) confirms the full
> divergence `∇·(δ∇(φ c_sat))` — the consistent form is correct. The original `DpDu` form is
> exactly what the Numerical Implementation doc §1.3.1.1/[δ^φ] prescribes, but that document's
> own appendix derivation (§2.14.3) starts from the strong first-order form and silently drops
> `φ∇·(δ∇c_sat)` relative to its eq. 515 — an erratum in the document, not a deliberate scheme.
> **Retracted:** the thermal-side vapor/liquid conduction terms (`Element2D.cxx:496-539`) do
> NOT need the same correction — per Theoretical eq. 424 they are advective enthalpy-transport
> terms (`C_l g_l·∇T`, `C_v g_v·∇T`), for which the strong-form `DpDu` treatment is the correct
> discretization. Decision: keep the consistent form in the moisture element only.

### D2. RH primary variable saturates (structural)
`φ ∈ [0,1]` cannot represent moisture above `w(φ=1)`. Near saturation `ξ = dw/dφ` blows up
(CottaerSandstone: `w` 63→180 kg/m³ over `φ` 0.99→1.0, `TestMaterials.hxx`), so `M` dominates
and the system stiffens/ill-conditions exactly in the high-moisture regime under test.

### D3. Convergence accepts non-solutions
`transientTimestep` returns `converged` from a solution-*change* metric only; the residual is
computed (`Domain.cxx:455`) but never gates acceptance. Oscillation midpoint-averaging
(`Domain.cxx:145-173`) and all-DOF clamping (`MoistureDomain.cxx:197-228`) both report
`converged = true`, and the middle level then accepts the first Δt that "converges"
(`Domain.cxx:299-317`) → dt-/parameter-dependent answers.

### D4. Single-pass staggered coupling
No per-timestep coupling iteration (`MultiDomain.cxx:63-136`) → first-order splitting error
that scales with the (adaptive) sub-step size.

### D5. Known-wrong phase-change term
Disabled with an author note that it is incorrect and "solver is not producing correct
results" (`Element2D.cxx:458`); must be reconciled with the thermal latent-heat coupling.

### D6. Non-conservative capacity → moisture-mass drift (found via the 1D reference)
The engine discretizes the storage term as `ξ(φ) ∂φ/∂t` with `ξ = dw/dφ` evaluated at the
current iterate (`Cap(dw/dφ)`, `Element2D.cxx:594`; lumped mass, `Elements2D.cxx:115`). This
**non-conservative** form does not conserve total moisture `∫w(φ)` when `ξ` varies steeply
over a step — precisely the near-saturation regime (Cottaer `w`: 63→180 kg/m³ over `φ`
0.99→1.0). The clean-room 1D reference reproduces this: with zero-flux boundaries (no moisture
may enter or leave), the naive `ξ φ_t` form **leaked ~6.2% of total moisture** over a 20-step
run, while the mass-conservative form (secant capacity `(w−w_old)/(φ−φ_old)` with lumped mass,
the Celia et al. 1990 "modified Picard" scheme) conserved mass to machine precision. Spurious
creation/destruction of moisture directly produces non-physical, parameter-sensitive
results in high-moisture cases.

**Fix direction:** switch the moisture capacity to the mass-conservative (Celia) form. The
`hygrothermfem_python` reference implements and verifies both forms and is the oracle for
this change.

---

## 6b. Fix progress (branch `high-moisture-solver-audit`)

| Defect | Status | Notes |
|---|---|---|
| **D6** capacity | **Implemented + verified** | `SorptionSecantCapacity` replaces the tangent `dw/dphi` in the moisture element. Engine closed-strip drift 6% → **1.5e-4**. Thermal green. New test: `tst/units/transient/MoistureMassConservation.unit.cxx`. |
| **D1** vapor weak form | **Implemented (correct); convergence caveat** | `QLEDpDuConsistentIntegrator2D` gives the consistent by-parts form of the temperature-gradient vapour term; used in the moisture element only, thermal assembly byte-identical (thermal green). The term is non-symmetric, so very tight tolerances converge harder: `LowHumidity + TightTolerance` now needs the D3 convergence work to avoid a non-convergence throw in the diagnostic sweep (regular suite has no throws). Correctness is best shown against the 1D reference under an imposed ∇T (todo). |
| **D3** convergence | **Implemented — relaxation-independent; tolerance/iteration + coupling residual** | Moisture converges on the **reduction of its free-DOF residual** to `1e-6 · ‖r₀‖` (relative to each problem's own initial residual, so it auto-scales — fixing the earlier fixed-threshold dead-end), with **best-effort acceptance** of the lowest-residual iterate instead of throwing, the change-metric oscillation exit disabled, and **adaptive under-relaxation** (halve damping when the residual grows; moisture ignores the user relaxation entirely and controls it itself). Moisture-only via `useResidualConvergence()`/`constrainedDofs()`; thermal byte-identical (green). **Sweep 14→9/42, 0 throws; LowHumidity/MediumHumidity fully consistent and HighHumidity consistent across all relaxation settings** — the relaxation-driven parameter-sensitivity is resolved. Remaining variability is in the tolerance / max-iteration knobs (best-effort accepts a less-converged iterate under `FewIterations`) and the staggered coupling for the extreme-ΔT `ExtremeHumidity` case → **D4**. The synthetic closed T-gradient diagnostic still leaks ~1% (iteration stalls there), so that residual-quality corner is open. |
| **D4** coupling | Identified as remaining lever | The residual `ExtremeHumidity` (80→20 °C) variability with tolerance/iteration settings traces to the single-pass staggered coupling (`MultiDomain::transient`) interacting with per-domain convergence. A per-timestep coupling (Picard) iteration would close it. |
| **D2** RH saturation | Smaller than thought — see diagnosis | Under D3's free-DOF residual convergence, `NearSaturation` went to **zero** violations across all settings — so it was a *convergence* problem (D3), not a representation problem. D2's real footprint is just the combined high-ΔT + near-saturation `ExtremeHumidity` case. |
| **D6-refine** capacity | **Implemented + verified** | Assemble the moisture capacity as the exact nodal-lumped secant diagonal `diag(vᵢ·ξᵢ)` (`m_LumpCapacityNodally`, moisture-only; thermal mass unchanged, thermal green). Closed-strip (uniform-T) conservation improved **1.5e-4 → 2.4e-15 (machine precision)**. Result *isolates the remaining leak*: the imposed-T-gradient case still leaks ~1% — but since the capacity is now provably exact, that leak is **not capacity; it is a residual imbalance from incomplete convergence → D3.** The reference, iterating to a tight residual, conserves that same case to 1e-12. |
| **D5** phase change | Deferred | Above-freezing cases unaffected. |

**Solver status now (D6+D1 on branch):** all pure-thermal tests green; **0 exceptions in the
regular suite**; 29 moisture/coupled golden tests are expected-red (clean value mismatches —
they encode pre-fix behaviour). The diagnostic sweep has one `TightTolerance` non-convergence
throw attributable to D1's non-symmetric term (awaiting the D3 redesign).

**Golden-test re-baseline (pending review):** the moisture goldens must be regenerated from the
fixed engine, but deliberately **not blind-snapshotted here** — the near-saturation values are
still formulation-limited (D2) and the D1 convergence caveat is open, so the new expected
values should be reviewed (and the isothermal moisture-only cases cross-checked against the 1D
reference, where D6 conservation + D1-inert make them the most trustworthy) before they are
frozen as golden.

## 6c. Continuation guide — pick up here

Written as a handoff so a later session can resume without re-deriving the above. Everything is
on branch `high-moisture-solver-audit` (local, **not pushed**). Companion: the clean-room 1D
reference at `~/Programming/hygrothermfem_python` (its own repo) — the **source of truth for
correctness** (the sweep is only a smoke test).

### Reference-document cross-check — DONE (2026-07-10)
The THERMM documents (Theoretical Model rev. 03/20/2020, Numerical Implementation rev.
5/28/2019, in `D:\tmp\htf numerical\`) were cross-checked against §2–§4 and the code:
(a) **Confirmed** — the vapour driving potential is `∇·(δ ∇(φ c_sat))` (Theoretical eq. 36,
Numerical eq. 515); `p_sat`, `R_v = 461.4`, `D_φ = ξ·D_l` and the moisture BC all match the
code exactly. (b) **D1 resolved as an erratum in the Numerical doc** — its `[δ^φ]` matrix has
the test function undifferentiated (the original `DpDu`), but its appendix derivation drops
`φ∇·(δ∇c_sat)` relative to its own eq. 515; the Theoretical PDE supports the consistent
by-parts form. Decision: keep `DpDuConsistent`. Thermal-side `DpDu` terms are advective
(Theoretical eq. 424) and correct as written. (c) **Nothing load-bearing missing** — air
advection (no air-pressure solver, TODO in code), gravity liquid term (doc drops it too),
rain `D_φs/D_φr` switch, sources, hysteresis, freezing (D5) are documented-but-deferred.
**New finding:** the docs prescribe `E(T) = (22.2 + 0.14·T_C)·1e-6` (Hagentoft); code used
constant `2.5e-5` (the 20 °C value). Decision: adopt E(T) before the golden re-baseline.
D6 secant capacity deviates from Numerical §1.3.1.5 (tangent) but realizes Theoretical
eq. 16 discretely — kept. Full errata/decision log:
`D:\Documents\HygroThermFEM Solver\moisture_solver_cross_check.typ`.

### What is done and committed (all thermal-green, moisture-only where noted)
- **D1** — consistent weak form of the temperature-gradient vapour term
  (`QLEDpDuConsistentIntegrator2D`, used only by `ElementMoistureLinear2D`). Oracle-verified:
  isothermal exact, gradient drives moisture hot→cold like the reference.
- **D6 + D6-refine** — mass-conservative capacity: secant `SorptionSecantCapacity` + exact
  nodal lumping (`m_LumpCapacityNodally`). Closed-strip conservation to **2.4e-15**.
- **D3** — settings-independent moisture convergence: free-DOF residual-*reduction*
  (`1e-6·‖r₀‖`) + best-effort acceptance of the lowest-residual iterate + adaptive
  under-relaxation (moisture ignores the user relaxation). Gated by `useResidualConvergence()`
  so thermal is byte-identical. Achieved relaxation-independence.
- **Sweep** — refactored to a mesh-refinement axis (3/6/12/24 elements) × tolerance/iteration
  settings; relaxation axis removed (solver-controlled now). Diagnostic only, **not in CI**.

### Test state (the standard gtest suite is what CI runs)
- **149 pass, 29 fail. All 29 failures are moisture/coupled golden tests; zero pure-thermal.**
  They hardcode pre-fix numbers, so the correctness fixes shifted them. Full list is the
  `Moisture*`, `MultiDomain*`, `SteadyState_2D_Exclude*`, `TestModelWithFrameCavity3`,
  `MoistureBC_2D_*` tests. **Next step to green CI: re-baseline the *water-content* golden
  arrays** (leave temperature arrays) from the fixed engine. Recommended: cross-check the
  isothermal moisture-only ones against the oracle first (most trustworthy); the near-saturation
  ones are still formulation-limited (D2), so eyeball them.

### Dead ends — do NOT repeat these
- **Subdivide-on-clamp / report-actual-convergence + subdivide**: makes non-converging cases
  subdivide to dt/1000 on *every* coupling sub-step and grind through thousands of tiny steps →
  runaway (minutes on tiny beams). Best-effort "always accept" exists precisely to avoid this.
- **Fixed absolute residual gate** (`‖r‖/‖B‖ ≤ tol`): no single threshold generalizes
  (near-saturation wants 1e-6, capacity-dominated LowHumidity can't reach it → throws). The
  *reduction* form (`‖r‖ ≤ 1e-6·‖r₀‖`) is what works.
- **Monotonic-only damping** (never recover toward 1.0): converges too slowly; 13/42 vs 9/42.
  Keep the halve-on-growth / ×1.5-recovery form.

### Known real defect still open: first-step overshoot
Coupled scenarios with a temperature jump drive humidity **negative → clamped to 0 (w=0) at
step 1 only** (40 such cases in the mesh sweep, all step 1). It's an initial-condition/boundary
shock the big first backward-Euler step doesn't resolve, and best-effort acceptance no longer
subdivides it. Needs a *lightweight* first-step resolver (e.g. subdivide only the first step, or
only when a DOF newly hits a bound, with a small bounded number of subdivisions) — the naive
full-subdivision approach is the runaway dead-end above. NOT caused by the multi-material
water-content calc (`Node2D::calcWaterContent` normalizes by total weight — verified correct;
single-material beams show the same overshoot).

### Known real defect still open: start-up thermal ringing (latent-coefficient feedback)
Found 2026-07-13 on the THERM "Ireggular_1" model (Stucco/Laminated panel/Fiberglass, saturated
start phi=0.95, constant time-series BCs 10C/0.1/7 vs 40C/0.3/10, dt=360): ~48 of 143 nodes
flip-flop in TEMPERATURE from step 2, max amplitude 9.3 K, damping to zero by step ~13.
Humidity is quiet (2 tiny flips) — thermal-side only. Flipping nodes cluster in FIBERGLASS
along the warm interior boundary. Signature is step-count-anchored, not physical-time-anchored
(same ~10-step ringing at dt=360 and dt=36; smooth at dt=3.6), i.e. algorithmic.
Hypothesis: temperature feedback of the (now-implicit) latent conductance
k_lat = h_lg * delta * phi * c_sat'(T) — near 40 C it is ~4x fiberglass dry k, and the
coefficient lags the solve between steps/iterations, giving damped over/undershoot ringing on
large warm-up increments. Candidate remedies: evaluate the latent DDu coefficient at the
current Picard iterate (verify it is not previous-timestep-lagged), under-relax the coefficient
update, or first-step subdivision (shares machinery with resolveShockStep).
Mesh refinement does NOT help (verified in THERM 2026-07-13) — consistent with the
time-stepping-anchored signature: the cause is temporal (coefficient lag), not spatial
resolution, so only dt/iteration-side remedies are candidates.
DRIVER CONFIRMED (THERM GUI bisection 2026-07-13): excluding Heat of Evaporation removes the
ringing entirely. Same latent term as the (fixed) explicit-assembly sawtooth, different
mechanism: the now-implicit T-part's COEFFICIENT k_lat(T, phi) is re-evaluated between
steps/iterations and lags the solve.
ENGINE REPRODUCTION: tst/units/therm_samples/ThermSample_Irregular1.unit.cxx
(DISABLED_StartupRinging; mesh embedded in Irregular1Mesh.hxx, generated from the THMZ by
D:/tmp/gen_irregular1_mesh.py). Reproduces THERM exactly: maxAmp 9.3365 K vs 9.3358 from the
THMZ analysis, last flip step 13 in both. Runs in ~2 s. Promotion criterion when fixed:
enable the test (asserts maxAmp < 0.5 K and no flips beyond step 3).
Detection tooling (works on any THERM transient THMZ, no GUI needed):
`D:\tmp\thmz_flip_detector.py <thmz>\transient results` (per-step flip counts + amplitudes)
and `D:\tmp\thmz_node_locator.py` (flip-node coordinates/materials). The full FE mesh is in
the THMZ's `transient results/Geometry.xml`, so the case can be rebuilt as an engine unit test
(arbitrary quads via createNode/createElement) for mesh/dt refinement without THERM.

### Remaining levers (priority order) — UPDATED 2026-07-10 evening
1. **Re-baseline moisture goldens** (33 red, all value-shift) → green CI → push.
2. ~~First-step overshoot~~ **RESOLVED**: the catastrophic step-1 collapse (φ→0 clamp) no
   longer reproduces after the conservation fixes — an instrumented sweep proved moisture
   solutions never hit bounds there anymore. Added as insurance anyway: a bounded
   new-bound-hit retry (`IDomain::resolveShockStep`, moisture-only via
   `retryStepOnNewBoundHit()`): when a converged step lands DOFs on a physical bound they
   did not start at, the step is redone ONCE as 8 fixed substeps (no recursion — not the
   subdivide-on-clamp runaway). Fires exactly once in the whole suite (a hard-drying
   golden case). The sweep's remaining "violations" were a NAIVE INVARIANT:
   `expectNonDecreasing(water)` applied to non-isothermal scenarios, where condensation
   on a cold surface followed by drying as it warms is legitimate physics — now gated to
   isothermal scenarios only. **Sweep after all fixes: 0 / 120 combinations with issues,
   0 throws** — full settings- and mesh-independence achieved.
3. **D4 coupling** — no longer visibly needed (NearSaturation/ExtremeHumidity variability
   gone from the sweep); revisit only if golden re-baseline cross-checks disagree.
4. **D2** — same status: only if evidence appears.

### How to verify (correctness, not the smoke test)
- Oracle self-checks: `cd ~/Programming/hygrothermfem_python && uv run pytest` (MMS order +
  machine-precision conservation).
- Engine-vs-oracle D1 isolation: build `HygroThermFEM_tests`, run
  `--gtest_filter=VaporTemperatureDiffusion.*` (writes `/tmp/htf_d1_*.csv`), then
  `uv run python examples/compare_d1.py`.
- Mass conservation: `--gtest_filter=MoistureMassConservation.*` (drift must be ~1e-15).

### Key files
- Equations/assembly: `src/elements/Element2D.cxx` (`ElementMoistureLinear2D`,
  `QLEDpDuConsistentIntegrator2D`, `nodalLumpedCapacity`), `src/math/Functions.cxx`
  (`SorptionSecantCapacity`, `SaturationFunction`).
- Solver/convergence (D3): `src/domain/Domain.cxx` (`transientTimestep`, `performNRIteration`,
  `evaluateConvergence`, `transient`), `src/domain/MoistureDomain.cxx`
  (`useResidualConvergence`, `constrainedDofs`).
- Coupling (D4): `src/domain/MultiDomain.cxx` (`transient`).
- Diagnostics/tests: `tst/units/transient/MoistureMassConservation.unit.cxx`,
  `VaporTemperatureDiffusion.unit.cxx`; `tst/sweep/SolverParameterSweep.cxx`.

### Session log 2026-07-10 — oracle rebuilt, E(T) adopted, engine leak root-caused
- **Oracle rebuilt** at `D:\Programming\Python\hygrothermfem_python` (uv project, typed,
  scipy): both D1 forms + both capacity forms selectable; banded O(n) assembly/solve
  (2001 nodes / 24 steps in 0.05 s); 11 self-checks in ~1 s incl. MMS order 2, analytic
  steady state `φ = C/c_sat` exact, independent scipy-BDF cross-integration
  (`solve_reference`), machine-precision conservation; pytest-timeout 60 s.
- **Oracle findings on D1 (strong form):** (a) NO driving force on uniform φ under ∇T
  (trial row sums annihilate constant φ — physically wrong); (b) 6.5e-3 closed-strip
  drift on nonuniform φ (consistent: 5.8e-16); (c) O(1) steady-state bias ~0.26 in φ
  (consistent: exact to machine).
- **E(T) adopted** (`vaporDiffusionCoefficientAtTemperature`, `VaporPermeability` in
  Functions.cxx; three Element2D.cxx sites). Test suite: 148 pass / 30 fail — all
  failures the known moisture/coupled goldens (one more than before E(T), as expected).
- **Engine-vs-oracle:** isothermal D1 diagnostic agrees EXACTLY (diff 0.0). Gradient
  case: engine *creates* 1.05% moisture monotonically and over-transports (up to 0.099
  in φ). Root cause (proved in the oracle, which had the identical signature): Picard
  residual stalls in a period-2 cycle because the liquid coefficient uses a pointwise
  tangent of the tabulated sorption curve — discontinuous at breakpoints; damping cannot
  break the cycle (fixed ω=0.5 still cycles). Oracle fix that removed it entirely:
  element-secant slope `(w_j−w_i)/(φ_j−φ_i)` as the liquid coefficient — exact for the
  element flux `−D_l ∇w` and continuous in nodal humidities.
- **RESOLVED (same session): the 1% leak had two causes, both fixed.** A vapor-only
  variant of the diagnostic (liquid off, T frozen ⇒ linear steps ⇒ iteration quality ruled
  out) showed the leak unchanged ⇒ operator-level, not Picard. (1) **Engine defect:**
  `IQLEIntegrator2D::integrate` applies coefficients as the pairwise nodal average
  `½(k_i+k_j)`, whose column sums do NOT vanish when the coefficient varies inside an
  element (c_sat(T) under ∇T) — a distributed spurious moisture source. Fix:
  `integrateInterpolated` (coefficient interpolated at Gauss points; column sums
  telescope to zero for ANY varying coefficient), enabled moisture-only via
  `m_InterpolateCoefficientsAtGaussPoints`; thermal byte-identical. Drift 1.05e-2 → 1.3e-3.
  (2) **Test-harness defect:** the diagnostics call `moisture().transient()` directly and
  never rolled the nodes' previous-timestep humidity forward (that happens in
  `MultiDomain::transient`), so the secant capacity measured against the INITIAL state —
  production unaffected. Fixed in both diagnostic tests. Drift 1.3e-3 → 2.3e-8; uniform-T
  conservation 2.6e-14. **Engine now matches the oracle to 2.2e-5** on the gradient
  diagnostic (full physics AND vapor-only). The earlier pairwise-secant-liquid hypothesis
  is retracted for the engine (oracle-only lesson; engine iteration converges fine).
- **Liquid-curve discontinuity RESOLVED (user decision 2026-07-10):** the log/geometric
  interpolation with a (0,0) anchor made D_l ≡ 0 across the whole first segment and then
  JUMP 0 → 1e-8 at w = 27 (a step discontinuity exactly where near-saturation cases
  iterate; the user confirms log-interp was an iteration-stability guess, not physics).
  Fix: `LiquidTransportationCurve` floors exact-zero table values at a physically
  negligible 1e-15 m²/s (`floorZeroValues`, Functions.cxx), so D_l decays continuously
  and steeply toward the anchor (Cottaer: ~1e-13 at w = 10) instead of the hard 0 + cliff.
  Oracle mirrors it (`MIN_LIQUID_DIFFUSIVITY`). Engine↔oracle agreement after: 2.2e-5.
  **Consequences:** (1) three more goldens shifted red — `SteadyState_2D_1`,
  `SolverSettingsInjection`, `MoistureBC_2D_2` (thermal capillary-conduction term now
  nonzero in the previously dead 0<w<27 zone; shifts ~2.5e-5 °C) → suite 146 pass /
  33 red, all value-shift goldens, re-baseline set grows to 33. (2) Sweep: NearSaturation
  and MediumHumidity violations went to ZERO (the w=27 cliff was hurting them — floor
  helps exactly where predicted), but step-1-overshoot-prone scenarios got noisier
  (51→65 combos with issues, still 0 throws; all new violations are step-1 water-content
  drops on coarse meshes — the KNOWN first-step overshoot expressing differently because
  transiently-dry iterates now see tiny nonzero D_l). ⇒ Recommended order change: fix the
  **first-step overshoot BEFORE re-baselining goldens**, so goldens are frozen once.
- **Suite after fixes: 149 pass / 30 fail — all 30 are the moisture/coupled goldens**
  (zero thermal), ready for re-baseline. Full write-up: typst cross-check doc.

## 7. What to verify next (Python 1D reference)

- **Scenario 1 — isothermal, moisture-only:** term (B) is inert, so C++ and a correct 1D
  solver should agree; isolates D2/D3.
- **Scenario 2 — imposed `∇T`:** term (B) is active; expect C++ to diverge from a consistent
  reference — this is the direct test of D1.
- **Scenario 3 — fully coupled:** exposes D4.
- **Oracle self-checks:** manufactured-solution convergence order and step-wise mass
  conservation (both independent of any external tool).
