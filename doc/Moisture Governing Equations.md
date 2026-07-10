# Moisture governing equations — reconstructed from the assembly code

> **Status: reverse-engineered, not authoritative.** There is no original
> differential-equation reference in the repo for the moisture model. This document
> reconstructs the continuous PDE, its weak form, and the discrete system **directly from the
> element-assembly code** so the equations can be audited. Where the code and a consistent
> finite-element derivation disagree, that is called out explicitly in
> [§6 Defects](#6-defects-found). Cross-check against the original ORNL/LBNL formulation when
> it surfaces.

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
([Phase B, scenario 2](#)). The analogous thermal-side vapor/liquid conduction terms
(`Element2D.cxx:496-539`) use the same `DpDu` pattern and likely need the same correction.

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
| **D3** convergence | **Attempted twice, reverted — hard; direction found** | Goal: make the accepted solution independent of the relaxation/tolerance/max-iteration settings. (a) A naive relative-residual gate fights the clamping (clamped DOFs carry an irreducible residual → false non-convergence → throws). (b) Fixed with a **free-DOF residual criterion** (residual over non-bound DOFs only, moisture-only via a `useResidualConvergence()` virtual so thermal stays byte-identical). This *worked* — at gate 1e-6 `NearSaturation` became perfectly settings-consistent (0 violations across all 7 settings). **But no single fixed threshold generalizes**: 1e-6 is unreachable for the capacity-dominated easy `LowHumidity` case (throws), while 1e-5 lets `LowHumidity` converge but is too loose for `NearSaturation`. The residual scale is problem-dependent. Correct next step: a **problem-scaled residual** (e.g. Jacobi/diagonal-scaled, or normalized by the mass/stiffness split) or **best-effort acceptance** (return the lowest-residual iterate at the finest subdivision instead of throwing). Reverted to keep the tree clean; the `useResidualConvergence`/`constrainedDofs` scaffolding is the intended shape. |
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

## 7. What to verify next (Python 1D reference)

- **Scenario 1 — isothermal, moisture-only:** term (B) is inert, so C++ and a correct 1D
  solver should agree; isolates D2/D3.
- **Scenario 2 — imposed `∇T`:** term (B) is active; expect C++ to diverge from a consistent
  reference — this is the direct test of D1.
- **Scenario 3 — fully coupled:** exposes D4.
- **Oracle self-checks:** manufactured-solution convergence order and step-wise mass
  conservation (both independent of any external tool).
