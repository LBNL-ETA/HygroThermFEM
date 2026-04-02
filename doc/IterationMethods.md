# Iteration Methods in HygroThermFEM

This document describes the numerical iteration strategies used in HygroThermFEM for solving coupled hygrothermal (heat and moisture) problems.

## 1. Problem Overview

HygroThermFEM solves two coupled partial differential equations over a 2D finite element domain:

- **Thermal equation**: governs temperature distribution, with material properties (thermal conductivity) that may depend on moisture content.
- **Moisture equation**: governs relative humidity distribution, with transport coefficients that depend on temperature.

Because material properties depend on the solution itself, both equations are in general **nonlinear**. The coupling between the two domains adds a further layer of nonlinearity: temperature affects moisture transport, and moisture affects thermal conductivity.

The solver addresses this through a hierarchy of three levels:

```
Outer level:   Adaptive time-accumulation with cross-coupling  (MultiDomain::transient)
Middle level:  Adaptive timestep probing                       (IDomain::transient)
Inner level:   Newton-Raphson iteration                        (IDomain::transientTimestep)
```

## 2. Newton-Raphson Iteration (Inner Level)

Each individual domain (thermal or moisture) solves its nonlinear system using a **modified Newton-Raphson (NR) method with relaxation**.

### 2.1 Algorithm

Given the discretized FEM system `A * U = B`, where both `A` (stiffness + mass + boundary matrices) and `B` (load vector) depend on the current solution `U`:

```
Input:  U_0 = initial guess (previous timestep solution)
        omega = relaxation parameter (default 1.0)
        tol = convergence tolerance (default 1e-5)

1.  Assemble A(U_0), B(U_0)
2.  For k = 0, 1, 2, ...
      a.  Compute residual:        r = B - A * U_k
      b.  Solve for correction:    A * dU = r
      c.  Limit increment:         allClamped = clamp dU per DOF (see Section 4.1)
      d.  Apply relaxed update:    U_{k+1} = U_k + omega * dU
      e.  Post-process:            enforce physical bounds on U_{k+1}
      f.  Update node properties with U_{k+1}
      g.  Reassemble:              A(U_{k+1}), B(U_{k+1})
      h.  Check convergence:       metric = |norm(U_{k+1}) - norm(U_k)| / norm(U_{k+1})
      i.  If metric <= tol: converged, exit
      i2. If allClamped: converged at physical bounds, exit (see Section 4.1)
      j.  If oscillation detected: apply midpoint averaging (see Section 4.2), exit
      k.  If k > MaxIterations: not converged, exit
```

The convergence criterion is based on the **relative change in the L2 norm** of the solution between consecutive iterations. This provides a global measure of whether the solution has stabilized.

### 2.2 Linear Problems

When both boundary conditions and element properties are linear (no state-dependent materials), the system `A * U = B` is solved in a single direct step without iteration.

## 3. Adaptive Time-Accumulation with Cross-Coupling (Outer Level)

The coupled hygrothermal problem is solved using a **partitioned scheme with adaptive time-stepping at the multi-domain level**. Rather than assembling a monolithic system, each domain is solved independently and then cross-domain data is exchanged after every successful sub-step.

### 3.1 Design Rationale

Earlier versions of the solver used a staggered coupling loop where each domain internally subdivided the timestep as needed, grinding through potentially thousands of tiny sub-steps with the cross-domain field (temperature or humidity) frozen at its initial value. For extreme conditions (high humidity at elevated temperatures), this led to:

- **Excessive computation**: the moisture domain might require 40,000+ internal sub-steps to cover a single 3600s timestep
- **Reduced accuracy**: moisture transport coefficients depend on temperature, but that temperature was frozen during all internal sub-steps
- **Poor coupling convergence**: the staggered outer loop had to repeat the full internal solve to capture cross-domain effects

The current approach moves the time-accumulation loop out of the individual domains and into `MultiDomain::transient`. Each domain now returns after its **first successful solve** (at whatever sub-timestep it could converge), and cross-coupling data is exchanged after every sub-step. This ensures both domains always see each other's latest solution.

### 3.2 Transient Coupling Algorithm

```
Input:  T_0 = initial temperature, H_0 = initial humidity, dt = target timestep

1.  totalTime = 0
2.  While totalTime < dt:
      a.  remainingTime = dt - totalTime
      b.  effectiveDt = remainingTime

      c.  Solve moisture domain at effectiveDt:
            - Moisture domain probes dt, subdivides if NR fails,
              returns first successful solution with actual dt used
            - effectiveDt = min(effectiveDt, moisture_dt)

      d.  Solve thermal domain at effectiveDt:
            - Thermal domain probes the (possibly reduced) dt
            - If thermal needs even smaller dt:
                effectiveDt = thermal_dt
                Re-solve moisture at effectiveDt

      e.  Accept both solutions, advance time:
            - T_current = T_new, H_current = H_new
            - Update nodes with both fields (cross-coupling exchange)
            - totalTime += effectiveDt

3.  Final update: advance "previous timestep" values in nodes
```

### 3.3 Key Properties

**Adaptive synchronization**: both domains always advance by the same sub-timestep. When the moisture domain (typically the stiffer one) needs a smaller dt, the thermal domain is solved at the same reduced dt. If the thermal domain requires an even smaller step (rare in practice), moisture is re-solved at that step.

**Frequent coupling exchange**: cross-domain data is refreshed after every successful sub-step, not after a full dt traversal. This dramatically improves convergence for strongly coupled problems: instead of solving moisture with a frozen temperature field, the temperature is updated at every sub-step.

**Performance improvement**: for the extreme test case (80C, RH=0.9999), the new approach reduced solve time from ~15.5 seconds (40,000 level-2 sub-steps with frozen coupling) to ~0.3 seconds (a few dozen sub-steps with live coupling) -- approximately 50x faster.

### 3.4 Steady-State Coupling

The steady-state coupling follows a simpler pattern: alternately solve each domain and exchange data, repeating until both converge. No time-accumulation is needed.

## 4. Convergence Stabilization Techniques

Near-saturation humidity conditions (RH close to 1.0) present particular challenges for convergence. Several stabilization techniques are employed:

### 4.1 Per-DOF Increment Limiting

**Problem**: The Newton-Raphson correction `dU` can project the humidity solution outside the physical range [0, 1]. Post-process clamping (forcing values back into bounds) introduces an inconsistency: the system matrices were assembled assuming an unconstrained state, but the actual solution is clamped. This causes the solver to repeatedly overshoot and clamp, creating artificial oscillations.

**Solution**: Before applying the correction, each DOF's increment is scaled so the projected value `U + omega * dU` remains within [0, 1]:

```
For each DOF i:
    projected = U[i] + dU[i] * omega
    if projected > 1.0:  dU[i] = (1.0 - U[i]) / omega
    if projected < 0.0:  dU[i] = -U[i] / omega
```

This keeps the solution within physical bounds at every iteration while maintaining consistency between the solution state and the assembled system matrices.

**Early convergence via saturation detection**: `limitIncrement` returns a boolean flag indicating whether **all** DOF corrections were clamped to effectively zero (absolute value below `1e-12`). When `allClamped` is true, the solution is pinned at physical bounds (e.g., humidity = 0 or 1 at every node) and no correction can improve it further. The NR loop treats this as converged immediately, avoiding unnecessary iterations when the solution is physically saturated. The base class (`IDomain`) returns `false` by default, so domains without physical bounds (e.g., thermal) are unaffected.

### 4.2 Oscillation Detection and Midpoint Averaging

**Problem**: Piecewise-linear material data (e.g., sorption isotherms) can have steep kinks where the derivative changes abruptly. When the NR iteration crosses such a kink, the Jacobian (effectively the assembled `A` matrix) changes discontinuously, causing the solution to alternate between two states on consecutive iterations (a "2-cycle"). The convergence metric stabilizes at a constant nonzero value instead of decreasing toward zero.

**Detection**: After a minimum warmup period (4 iterations), the solver compares the convergence metric between consecutive iterations. If the relative change is less than 1%:

```
|metric_k - metric_{k-1}| / metric_{k-1} < 0.01
```

the solver concludes that a 2-cycle is present.

**Resolution**: The solution is set to the **midpoint of the last two iterates**:

```
U_final = 0.5 * (U_k + U_{k-1})
```

This places the solution at the center of the oscillation band, which lies at (or very near) the material property kink. The midpoint is the best physically meaningful estimate when the true solution sits exactly on a non-smooth point of the material curve.

### 4.3 Smooth Vapor Transfer Coefficient

**Problem**: The water vapor transfer coefficient `beta` originally had a hard step discontinuity at RH = 1.0:

```
beta = value   if RH <= 1.0
beta = 0       if RH > 1.0
```

During NR iteration, small overshoots past RH = 1.0 would zero out the transfer coefficient, then on the next iteration (after clamping back to 1.0) it would jump back to full value, compounding oscillations.

**Solution**: Replace the step function with a smooth linear ramp:

```
smoothFactor = clamp((1.0 + w - RH) / w, 0, 1)     where w = 0.01
beta = betaValue * smoothFactor
```

This produces: full transfer below RH = 1.0, a linear taper over [1.0, 1.01], and zero transfer above RH = 1.01. The narrow transition width ensures negligible impact on results while eliminating the Jacobian discontinuity.

## 5. Adaptive Timestep Probing (Middle Level)

When the NR iteration fails to converge at a requested timestep, the domain does not give up. Instead, it **subdivides** and tries progressively smaller timesteps until one succeeds, then **returns immediately** with the first successful solution and the actual dt used:

```
Input:  dt = requested timestep, maxLevels = 3, numSub = 10

1.  currentDt = dt, level = 0
2.  Loop:
      a.  Attempt NR solve for currentDt
      b.  If converged:
            - Return (solution, currentDt)   <-- actual dt may be < requested dt
      c.  If not converged:
            - Subdivide: currentDt = currentDt / numSub
            - level += 1
            - If level > maxLevels: throw error (solution failed)
```

The key difference from earlier versions: the domain **does not accumulate time internally**. It finds the largest dt at which the NR solver converges and returns that single-step result. The caller (`MultiDomain::transient`) is responsible for accumulating sub-steps until the target time is reached, exchanging cross-domain coupling data between each sub-step (see Section 3).

With default settings (3 levels, 10 subdivisions per level), the solver can probe timesteps as small as `dt / 1000` before reporting failure. Each subdivision level is reported to registered observers for progress monitoring.

## 6. Diagnostic Output

Both `IDomain` and `MultiDomain` support an optional diagnostic stream (`setDiagnosticStream(std::ostream*)`). When enabled, the solver writes CSV-formatted data for every NR iteration and sub-step:

**NR iteration columns**: `subtimestep_dt, nr_iter, metric, converge_tol, correction_norm, damping_factor, all_clamped, converged, reason, u0, u1, ...`

**Comment lines** (`# ...`) mark sub-step boundaries, domain switches (MOISTURE/THERMAL), and time-accumulation progress.

Usage in a test:
```cpp
std::ofstream diagFile("solver_diagnostics.csv");
multiDomain.setDiagnosticStream(&diagFile);
```

Load in Python with: `pd.read_csv("solver_diagnostics.csv", comment="#")`

## 7. Summary of Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Relaxation parameter | 1.0 | NR under-relaxation factor (omega) |
| Error tolerance | 1e-5 | Convergence threshold for relative norm change |
| Max iterations | 50 | Maximum NR iterations per timestep |
| Max division levels | 3 | Maximum timestep probing depth |
| Subdivisions per level | 10 | Factor by which dt is reduced per probing level |
| Oscillation check threshold | 0.01 | Relative metric change below which 2-cycle is detected |
| Min iterations for oscillation check | 4 | Warmup before oscillation detection activates |
| Clamp tolerance | 1e-12 | Absolute increment below which a DOF is considered clamped |
| Vapor transfer transition width | 0.01 | RH range over which beta tapers to zero |

## 8. Licensing and Prior Art

All numerical methods used in this solver are standard techniques from the public domain, with no proprietary or patent-encumbered algorithms:

- **Newton-Raphson iteration** — published by Isaac Newton (1685) and Joseph Raphson (1690). Standard textbook method for nonlinear systems.
- **Under-relaxation** — classical stabilization technique, widely described in numerical methods literature.
- **Partitioned (staggered) coupling** — standard approach in multi-physics FEM solvers since the 1970s.
- **Adaptive time-stepping with subdivision** — ubiquitous in ODE/PDE solvers; standard numerical practice.
- **Per-DOF increment limiting** — standard practice in constrained nonlinear solvers.
- **Midpoint averaging for oscillation resolution** — basic numerical stabilization, not a proprietary technique.
- **Caller-level time-accumulation with frequent coupling exchange** — an engineering design choice for organizing cross-domain data exchange; built from the standard methods above.

The specific combination, tuning parameters, and implementation details are original engineering work by the HygroThermFEM authors.
