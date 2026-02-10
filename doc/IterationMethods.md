# Iteration Methods in HygroThermFEM

This document describes the numerical iteration strategies used in HygroThermFEM for solving coupled hygrothermal (heat and moisture) problems.

## 1. Problem Overview

HygroThermFEM solves two coupled partial differential equations over a 2D finite element domain:

- **Thermal equation**: governs temperature distribution, with material properties (thermal conductivity) that may depend on moisture content.
- **Moisture equation**: governs relative humidity distribution, with transport coefficients that depend on temperature.

Because material properties depend on the solution itself, both equations are in general **nonlinear**. The coupling between the two domains adds a further layer of nonlinearity: temperature affects moisture transport, and moisture affects thermal conductivity.

The solver addresses this through a hierarchy of three nested iteration levels:

```
Outer level:   Staggered coupling loop      (MultiDomain)
Middle level:  Adaptive timestep subdivision (IDomain::transient)
Inner level:   Newton-Raphson iteration      (IDomain::transientTimestep)
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

## 3. Staggered Coupling Loop (Outer Level)

The coupled hygrothermal problem is solved using a **partitioned (staggered) scheme**. Rather than assembling and solving a single monolithic system for both temperature and humidity simultaneously, each domain is solved independently and the cross-domain data is exchanged iteratively.

### 3.1 Transient Coupling Algorithm

```
Input:  T_0 = initial temperature, H_0 = initial humidity

1.  Do
      a.  Inner loop (solve each domain until individually converged):
            - Solve moisture domain -> H_new
            - Solve thermal domain  -> T_new
            - Repeat until both domain errors < tolerance
              OR 2-cycle stall detected (see Section 3.2)

      b.  Cross-coupling exchange:
            - Update node temperatures with T_new
            - Re-solve moisture domain with updated temperatures -> H_new
            - Update node humidities with H_new
            - Re-solve thermal domain with updated humidities   -> T_new

      c.  Increment outer iteration counter

2.  While (T_error > tol OR H_error > tol) AND iterations < max
```

The inner loop uses `AND` for its convergence check (exits when either domain converges, since further uncoupled iteration is wasteful). The outer loop uses `OR` (continues until **both** domains have converged), ensuring the coupled solution is self-consistent.

### 3.2 Inner Loop 2-Cycle Stall Detection

In strongly coupled problems near material saturation (e.g., humidity close to 1.0 at elevated temperatures), piecewise-linear material data causes the inner loop's Gauss-Seidel iteration to oscillate between two states on consecutive iterations — a "2-cycle". The coupling error alternates high/low without decreasing, so neither domain converges individually and the `AND` exit condition is never satisfied.

Rather than using a fixed iteration cap (which would starve larger models that need more inner iterations to converge), the inner loop uses **adaptive 2-cycle stall detection**. After a warmup period of 12 iterations, the solver compares the current max coupling error with the error from 2 iterations ago (the same phase of any oscillation cycle):

```
maxErr = max(humidityError, temperatureError)

if innerIterCount >= 12
   AND |maxErr - maxErr_{k-2}| / maxErr_{k-2} < 0.01:
       exit inner loop
```

When the same-phase error has changed by less than 1%, the inner loop is stuck in a 2-cycle and further iterations are wasteful. The loop exits and lets the outer loop's cross-coupling exchange (step 1b) drive convergence.

This approach is adaptive:
- **Non-stalling problems**: errors decrease iteration to iteration, so the stall condition is never triggered and the loop runs to convergence via the `AND` exit condition — identical to having no cap at all.
- **2-cycling problems**: the stall is detected after the warmup period and the loop exits, avoiding infinite iteration while preserving enough inner iterations for the solution state to develop.

### 3.3 Steady-State Coupling

The steady-state coupling follows a simpler pattern: alternately solve each domain and exchange data, repeating until both converge.

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

## 5. Adaptive Timestep Subdivision (Middle Level)

When the NR iteration fails to converge within the maximum number of iterations for a given timestep, the solver does not immediately give up. Instead, it **subdivides the timestep** and attempts to reach the target time through multiple smaller steps:

```
Input:  dt = requested timestep, maxLevels = 3, numSub = 10

1.  totalTime = 0, currentDt = dt, level = 0
2.  While totalTime < dt:
      a.  Attempt NR solve for currentDt
      b.  If converged:
            - Advance: totalTime += currentDt
            - Update state variables
      c.  If not converged:
            - Subdivide: currentDt = currentDt / numSub
            - level += 1
            - If level > maxLevels: throw error (solution failed)
```

With default settings (3 levels, 10 subdivisions per level), the solver can attempt timesteps as small as `dt / 1000` before reporting failure. Each subdivision level is reported to registered observers for progress monitoring.

## 6. Summary of Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Relaxation parameter | 1.0 | NR under-relaxation factor (omega) |
| Error tolerance | 1e-5 | Convergence threshold for relative norm change |
| Max iterations | 50 | Maximum NR iterations per timestep; also outer coupling loop limit |
| Inner stall warmup | 12 | Minimum inner iterations before 2-cycle stall detection activates (see Section 3.2) |
| Inner stall threshold | 0.01 | Relative same-phase error change below which inner loop 2-cycle is detected |
| Max division levels | 3 | Maximum timestep subdivision depth |
| Subdivisions per level | 10 | Number of sub-timesteps per subdivision |
| Oscillation check threshold | 0.01 | Relative metric change below which 2-cycle is detected |
| Min iterations for oscillation check | 4 | Warmup before oscillation detection activates |
| Clamp tolerance | 1e-12 | Absolute increment below which a DOF is considered clamped |
| Vapor transfer transition width | 0.01 | RH range over which beta tapers to zero |
