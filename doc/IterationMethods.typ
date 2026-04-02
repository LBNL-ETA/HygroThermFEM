#set document(title: "Iteration Methods in HygroThermFEM", author: "LBNL / ORNL")
#set page(margin: 2.5cm, numbering: "1")
#set text(font: "New Computer Modern", size: 11pt)
#set heading(numbering: "1.1")
#set par(justify: true)

#show heading.where(level: 1): it => {
  pagebreak(weak: true)
  v(0.5em)
  text(size: 16pt, weight: "bold", it)
  v(0.4em)
}

#show heading.where(level: 2): it => {
  v(0.8em)
  text(size: 13pt, weight: "bold", it)
  v(0.3em)
}

#show heading.where(level: 3): it => {
  v(0.5em)
  text(size: 11pt, weight: "bold", it)
  v(0.2em)
}

#show raw.where(block: true): it => {
  set text(size: 9pt)
  block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%, it)
}

#align(center)[
  #text(size: 20pt, weight: "bold")[Iteration Methods in HygroThermFEM]
  #v(0.5em)
  #text(size: 12pt, fill: luma(100))[Lawrence Berkeley National Laboratory / Oak Ridge National Laboratory]
  #v(2em)
]

= Problem Overview

HygroThermFEM solves two coupled partial differential equations over a 2D finite element domain:

- *Thermal equation*: governs temperature distribution, with material properties (thermal conductivity) that may depend on moisture content.
- *Moisture equation*: governs relative humidity distribution, with transport coefficients that depend on temperature.

Because material properties depend on the solution itself, both equations are in general *nonlinear*. The coupling between the two domains adds a further layer of nonlinearity: temperature affects moisture transport, and moisture affects thermal conductivity.

The solver addresses this through a hierarchy of three levels:

#figure(
  table(
    columns: (auto, auto),
    align: left,
    stroke: 0.5pt,
    inset: 8pt,
    [*Level*], [*Responsibility*],
    [Outer: `MultiDomain::transient`], [Adaptive time-accumulation with cross-coupling],
    [Middle: `IDomain::transient`], [Adaptive timestep probing],
    [Inner: `IDomain::transientTimestep`], [Newton--Raphson iteration],
  ),
  caption: [Solver hierarchy],
)

= Newton--Raphson Iteration (Inner Level)

Each individual domain (thermal or moisture) solves its nonlinear system using a *modified Newton--Raphson (NR) method with relaxation*.

== Algorithm

Given the discretized FEM system $bold(A) dot bold(U) = bold(B)$, where both $bold(A)$ (stiffness + mass + boundary matrices) and $bold(B)$ (load vector) depend on the current solution $bold(U)$:

#block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%)[
  *Input:* $bold(U)_0$ = initial guess, $omega$ = relaxation parameter, $"tol"$ = convergence tolerance

  + Assemble $bold(A)(bold(U)_0)$, $bold(B)(bold(U)_0)$
  + For $k = 0, 1, 2, dots$
    + Compute residual: $bold(r) = bold(B) - bold(A) dot bold(U)_k$
    + Solve for correction: $bold(A) dot Delta bold(U) = bold(r)$
    + Limit increment: clamp $Delta bold(U)$ per DOF (see @sec-clamping)
    + Apply relaxed update: $bold(U)_(k+1) = bold(U)_k + omega dot Delta bold(U)$
    + Post-process: enforce physical bounds on $bold(U)_(k+1)$
    + Update node properties with $bold(U)_(k+1)$
    + Reassemble: $bold(A)(bold(U)_(k+1))$, $bold(B)(bold(U)_(k+1))$
    + Check convergence: $"metric" = |norm(bold(U)_(k+1)) - norm(bold(U)_k)| / norm(bold(U)_(k+1))$
    + If $"metric" <= "tol"$: converged, exit
    + If all DOFs clamped: converged at physical bounds, exit (see @sec-clamping)
    + If oscillation detected: apply midpoint averaging (see @sec-oscillation), exit
    + If $k >$ MaxIterations: not converged, exit
]

The convergence criterion is based on the *relative change in the $L^2$ norm* of the solution between consecutive iterations.

== Linear Problems

When both boundary conditions and element properties are linear (no state-dependent materials), the system $bold(A) dot bold(U) = bold(B)$ is solved in a single direct step without iteration.

= Adaptive Time-Accumulation with Cross-Coupling (Outer Level)

The coupled hygrothermal problem is solved using a *partitioned scheme with adaptive time-stepping at the multi-domain level*. Rather than assembling a monolithic system, each domain is solved independently and then cross-domain data is exchanged after every successful sub-step.

== Design Rationale

Earlier versions of the solver used a staggered coupling loop where each domain internally subdivided the timestep as needed, grinding through potentially thousands of tiny sub-steps with the cross-domain field (temperature or humidity) frozen at its initial value. For extreme conditions (high humidity at elevated temperatures), this led to:

- *Excessive computation*: the moisture domain might require 40,000+ internal sub-steps to cover a single 3600#h(0.15em)s timestep
- *Reduced accuracy*: moisture transport coefficients depend on temperature, but that temperature was frozen during all internal sub-steps
- *Poor coupling convergence*: the staggered outer loop had to repeat the full internal solve to capture cross-domain effects

The current approach moves the time-accumulation loop out of the individual domains and into `MultiDomain::transient`. Each domain now returns after its *first successful solve* (at whatever sub-timestep it could converge), and cross-coupling data is exchanged after every sub-step. This ensures both domains always see each other's latest solution.

== Transient Coupling Algorithm

#block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%)[
  *Input:* $bold(T)_0$ = initial temperature, $bold(H)_0$ = initial humidity, $Delta t$ = target timestep

  + $t_"total" = 0$
  + While $t_"total" < Delta t$:
    + $t_"remain" = Delta t - t_"total"$, #h(1em) $Delta t_"eff" = t_"remain"$
    + Solve moisture domain at $Delta t_"eff"$:
      - Domain probes $Delta t_"eff"$, subdivides if NR fails, returns first successful solution with actual $Delta t$ used
      - $Delta t_"eff" = min(Delta t_"eff", Delta t_"moisture")$
    + Solve thermal domain at $Delta t_"eff"$:
      - If thermal needs even smaller $Delta t$: set $Delta t_"eff" = Delta t_"thermal"$, re-solve moisture
    + Accept both solutions, advance time:
      - $bold(T)_"current" = bold(T)_"new"$, #h(1em) $bold(H)_"current" = bold(H)_"new"$
      - Update nodes with both fields (cross-coupling exchange)
      - $t_"total" += Delta t_"eff"$
  + Final update: advance "previous timestep" values in nodes
]

== Key Properties

*Adaptive synchronization.* Both domains always advance by the same sub-timestep. When the moisture domain (typically the stiffer one) needs a smaller $Delta t$, the thermal domain is solved at the same reduced $Delta t$. If the thermal domain requires an even smaller step (rare in practice), moisture is re-solved at that step.

*Frequent coupling exchange.* Cross-domain data is refreshed after every successful sub-step, not after a full $Delta t$ traversal. This dramatically improves convergence for strongly coupled problems: instead of solving moisture with a frozen temperature field, the temperature is updated at every sub-step.

*Performance improvement.* For the extreme test case (80~°C, RH=0.9999), the new approach reduced solve time from ~15.5~s (40,000 level-2 sub-steps with frozen coupling) to ~0.3~s (a few dozen sub-steps with live coupling) --- approximately 50× faster.

== Steady-State Coupling

The steady-state coupling follows a simpler pattern: alternately solve each domain and exchange data, repeating until both converge. No time-accumulation is needed.

= Convergence Stabilization Techniques

Near-saturation humidity conditions (RH close to 1.0) present particular challenges for convergence. Several stabilization techniques are employed.

== Per-DOF Increment Limiting <sec-clamping>

*Problem.* The Newton--Raphson correction $Delta bold(U)$ can project the humidity solution outside the physical range $[0, 1]$. Post-process clamping introduces an inconsistency: the system matrices were assembled assuming an unconstrained state, but the actual solution is clamped. This causes the solver to repeatedly overshoot and clamp, creating artificial oscillations.

*Solution.* Before applying the correction, each DOF's increment is scaled so the projected value $U_i + omega dot Delta U_i$ remains within $[0, 1]$:

$ "For each DOF" i: quad cases(
  Delta U_i = (1 - U_i) / omega & "if" U_i + Delta U_i dot omega > 1,
  Delta U_i = -U_i / omega & "if" U_i + Delta U_i dot omega < 0,
) $

*Early convergence via saturation detection.* `limitIncrement` returns a flag indicating whether *all* DOF corrections were clamped to effectively zero ($|Delta U_i| < 10^(-12)$). When true, the solution is pinned at physical bounds and no correction can improve it further. The NR loop treats this as converged immediately.

== Oscillation Detection and Midpoint Averaging <sec-oscillation>

*Problem.* Piecewise-linear material data (e.g., sorption isotherms) can have steep kinks where the derivative changes abruptly. When the NR iteration crosses such a kink, the Jacobian changes discontinuously, causing the solution to alternate between two states on consecutive iterations (a "2-cycle").

*Detection.* After a minimum warmup period (4 iterations), the solver compares the convergence metric between consecutive iterations. If the relative change is less than 1%:

$ |"metric"_k - "metric"_(k-1)| / "metric"_(k-1) < 0.01 $

the solver concludes that a 2-cycle is present.

*Resolution.* The solution is set to the midpoint of the last two iterates:

$ bold(U)_"final" = 1/2 (bold(U)_k + bold(U)_(k-1)) $

This places the solution at the center of the oscillation band, which lies at (or very near) the material property kink.

== Smooth Vapor Transfer Coefficient

*Problem.* The water vapor transfer coefficient $beta$ originally had a hard step discontinuity at RH = 1.0. During NR iteration, small overshoots past RH = 1.0 would zero out the transfer coefficient, then on the next iteration it would jump back to full value, compounding oscillations.

*Solution.* Replace the step function with a smooth linear ramp:

$ "smoothFactor" = "clamp"((1 + w - "RH") / w, space 0, space 1), quad w = 0.01 $
$ beta = beta_"value" dot "smoothFactor" $

This produces: full transfer below RH = 1.0, a linear taper over $[1.0, 1.01]$, and zero transfer above RH = 1.01.

= Adaptive Timestep Probing (Middle Level)

When the NR iteration fails to converge at a requested timestep, the domain *subdivides* and tries progressively smaller timesteps until one succeeds, then *returns immediately* with the first successful solution and the actual $Delta t$ used:

#block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%)[
  *Input:* $Delta t$ = requested timestep, maxLevels = 3, numSub = 10

  + $Delta t_"current" = Delta t$, #h(1em) level = 0
  + Loop:
    + Attempt NR solve for $Delta t_"current"$
    + If converged: return (solution, $Delta t_"current"$) #h(2em) #text(fill: luma(120))[$arrow.l$ actual $Delta t$ may be $<$ requested]
    + If not converged: $Delta t_"current" = Delta t_"current" / "numSub"$, #h(1em) level += 1
    + If level $>$ maxLevels: throw error
]

The domain *does not accumulate time internally*. It finds the largest $Delta t$ at which the NR solver converges and returns that single-step result. The caller (`MultiDomain::transient`) is responsible for accumulating sub-steps until the target time is reached, exchanging cross-domain coupling data between each sub-step (see Section 3).

With default settings (3 levels, 10 subdivisions per level), the solver can probe timesteps as small as $Delta t slash 1000$ before reporting failure.

= Diagnostic Output

Both `IDomain` and `MultiDomain` support an optional diagnostic stream via `setDiagnosticStream`. When enabled, the solver writes CSV-formatted data for every NR iteration and sub-step.

*NR iteration columns:* `subtimestep_dt`, `nr_iter`, `metric`, `converge_tol`, `correction_norm`, `damping_factor`, `all_clamped`, `converged`, `reason`, `u0`, `u1`, ...

*Comment lines* (`#`) mark sub-step boundaries, domain switches, and time-accumulation progress.

```cpp
std::ofstream diagFile("solver_diagnostics.csv");
multiDomain.setDiagnosticStream(&diagFile);
```

Load in Python: `pd.read_csv("solver_diagnostics.csv", comment="#")`

= Summary of Parameters

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, center, left),
    stroke: 0.5pt,
    inset: 8pt,
    [*Parameter*], [*Default*], [*Description*],
    [Relaxation parameter], [$1.0$], [NR under-relaxation factor ($omega$)],
    [Error tolerance], [$10^(-5)$], [Convergence threshold for relative norm change],
    [Max iterations], [$50$], [Maximum NR iterations per timestep],
    [Max division levels], [$3$], [Maximum timestep probing depth],
    [Subdivisions per level], [$10$], [Factor by which $Delta t$ is reduced per probing level],
    [Oscillation check threshold], [$0.01$], [Relative metric change below which 2-cycle is detected],
    [Min iterations for oscillation check], [$4$], [Warmup before oscillation detection activates],
    [Clamp tolerance], [$10^(-12)$], [Absolute increment below which a DOF is considered clamped],
    [Vapor transfer transition width], [$0.01$], [RH range over which $beta$ tapers to zero],
  ),
  caption: [Solver configuration parameters],
)
