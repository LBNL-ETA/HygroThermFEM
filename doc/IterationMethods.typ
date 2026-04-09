#import "@preview/cetz:0.3.4"

#set document(title: "Iteration Methods in HygroThermFEM", author: "Simon Vidanovic")
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

// ─────────────────────────────────────────────────────────
// Title
// ─────────────────────────────────────────────────────────
#align(center)[
  #text(size: 20pt, weight: "bold")[Iteration Methods in HygroThermFEM]
  #v(0.5em)
  #text(size: 12pt, fill: luma(100))[Lawrence Berkeley National Laboratory]
  #v(2em)
]

// ═════════════════════════════════════════════════════════
= Problem Overview
// ═════════════════════════════════════════════════════════

HygroThermFEM solves two coupled partial differential equations over a 2D finite element domain:

- *Thermal equation*: governs temperature distribution $T(x,y,t)$, with material properties (thermal conductivity, latent heat effects) that may depend on moisture content.
- *Moisture equation*: governs relative humidity distribution $phi(x,y,t)$, with transport coefficients (vapor diffusion, liquid transport) that depend on temperature.

Because material properties depend on the solution itself, both equations are in general *nonlinear*. The coupling between the two domains adds a further layer of nonlinearity: temperature affects moisture transport, and moisture affects thermal conductivity.

== Why Coupling Matters

Consider a concrete wall initially at $T = 20 °C$ and $phi = 0.5$. When exposed to cold, humid outdoor air:

+ The temperature drop increases relative humidity near the surface (saturation concentration $c_"sat"(T)$ decreases with temperature).
+ Higher humidity drives moisture transport inward, releasing latent heat of condensation.
+ The released heat modifies the temperature field, which in turn changes transport coefficients.

This circular dependency means neither equation can be solved independently --- the solver must iterate between them until a self-consistent solution is reached.

== Solver Architecture

The solver addresses the coupled nonlinear problem through a hierarchy of three levels:

#figure(
  table(
    columns: (auto, auto, auto),
    align: left,
    stroke: 0.5pt,
    inset: 8pt,
    [*Level*], [*Class*], [*Responsibility*],
    [Outer], [`MultiDomain::transient`], [Adaptive time-accumulation with cross-coupling],
    [Middle], [`IDomain::transient`], [Adaptive timestep probing],
    [Inner], [`IDomain::transientTimestep`], [NR loop orchestration],
    [Iteration], [`IDomain::performNRIteration`], [Single NR step: correction, line search, convergence],
  ),
  caption: [Solver hierarchy. Each level handles a different aspect of the nonlinear coupled problem.],
)

The iteration level performs a single NR step: computing the correction, running the backtracking line search, and checking convergence. The inner level orchestrates the NR loop, bundling state into an `NRLoopState` struct. The middle level handles convergence failure by trying smaller timesteps. The outer level handles the coupling between temperature and moisture.

#figure(
  cetz.canvas(length: 1cm, {
    import cetz.draw: *

    // Outer box: MultiDomain time accumulation
    rect((-0.5, -0.3), (12.5, 7.8), stroke: (paint: rgb("#2266aa"), thickness: 1.5pt), radius: 6pt)
    content((6, 7.4), text(size: 9pt, weight: "bold", fill: rgb("#2266aa"))[Outer: MultiDomain time accumulation])

    // Sub-step boxes
    for (idx, xoff) in ((0, 0), (1, 4.2), (2, 8.4)) {
      let label = if idx < 2 { "Sub-step " + str(idx + 1) } else { "..." }

      rect((xoff, 0), (xoff + 3.8, 6.8), stroke: (paint: luma(120), thickness: 0.8pt), radius: 4pt, fill: luma(248))

      content((xoff + 1.9, 6.4), text(size: 8pt, weight: "bold")[#label])

      if idx < 2 {
        // Moisture block
        rect((xoff + 0.3, 3.5), (xoff + 3.5, 6.0), stroke: 0.5pt, fill: rgb("#e8f0ff"), radius: 3pt)
        content((xoff + 1.9, 5.6), text(size: 7pt, weight: "bold", fill: rgb("#336699"))[Moisture domain])
        content((xoff + 1.9, 5.0), text(size: 6.5pt)[Probe $Delta t$])
        content((xoff + 1.9, 4.5), text(size: 6.5pt)[NR iterate])
        content((xoff + 1.9, 4.0), text(size: 6.5pt)[Return first success])

        // Thermal block
        rect((xoff + 0.3, 0.5), (xoff + 3.5, 3.0), stroke: 0.5pt, fill: rgb("#fff0e8"), radius: 3pt)
        content((xoff + 1.9, 2.6), text(size: 7pt, weight: "bold", fill: rgb("#996633"))[Thermal domain])
        content((xoff + 1.9, 2.0), text(size: 6.5pt)[Probe $Delta t$])
        content((xoff + 1.9, 1.5), text(size: 6.5pt)[NR iterate])
        content((xoff + 1.9, 1.0), text(size: 6.5pt)[Return first success])

        // Arrow: moisture -> thermal (coupling)
        line((xoff + 1.9, 3.5), (xoff + 1.9, 3.0),
             stroke: (paint: rgb("#44aa44"), thickness: 1pt),
             mark: (end: ">", fill: rgb("#44aa44")))
      } else {
        content((xoff + 1.9, 3.5), text(size: 14pt, fill: luma(150))[...])
      }
    }

    // Arrows between sub-steps (coupling exchange)
    for xoff in (3.8, 8.0) {
      line((xoff + 0.05, 4.5), (xoff + 0.35, 4.5),
           stroke: (paint: rgb("#44aa44"), thickness: 1.2pt),
           mark: (end: ">", fill: rgb("#44aa44")))
      line((xoff + 0.05, 2.0), (xoff + 0.35, 2.0),
           stroke: (paint: rgb("#44aa44"), thickness: 1.2pt),
           mark: (end: ">", fill: rgb("#44aa44")))
    }

    // Legend
    content((6, -0.7), text(size: 7pt, fill: rgb("#44aa44"))[Green arrows = cross-coupling data exchange ($T arrow.l.r phi$)])
  }),
  caption: [Solver flow. Each sub-step solves moisture then thermal, exchanging coupling data between domains and between sub-steps. The outer loop accumulates sub-steps until the target $Delta t$ is reached.],
) <fig-flow>

// ═════════════════════════════════════════════════════════
= Newton--Raphson Iteration (Inner Level)
// ═════════════════════════════════════════════════════════

Each individual domain (thermal or moisture) solves its nonlinear system using a *modified Newton--Raphson (NR) method with relaxation*.

== The Nonlinear FEM System

After spatial discretization with bilinear quadrilateral elements, the transient FEM equation for a single domain takes the form:

$ (bold(M) / (Delta t) + bold(K) + bold(H)) dot bold(U)^((n+1)) = bold(M) / (Delta t) dot bold(U)^((n)) + bold(R) $

where:
- $bold(M)$ is the lumped mass (capacitance) matrix --- stores energy or moisture
- $bold(K)$ is the conductance (stiffness) matrix --- diffusion/conduction
- $bold(H)$ is the boundary condition matrix (convective terms)
- $bold(R)$ is the load vector (boundary condition driving terms)
- $bold(U)^((n))$ is the solution at the current timestep
- $bold(U)^((n+1))$ is the solution at the next timestep

For notational convenience, we write this as $bold(A) dot bold(U) = bold(B)$, where both $bold(A)$ and $bold(B)$ depend on the current estimate of $bold(U)$ (through state-dependent material properties).

== Why Material Properties Create Nonlinearity

For the moisture equation, the key material properties are the *sorption isotherm* and the *liquid transport coefficient*. The sorption isotherm $w(phi)$ relates water content to relative humidity. For Cottaer Sandstone (used in the test suite), the sorption curve is:

#figure(
  table(
    columns: 9,
    align: center,
    stroke: 0.5pt,
    inset: 5pt,
    [$phi$], [0.0], [0.5], [0.8], [0.93], [0.95], [0.99], [0.999], [1.0],
    [$w$ (kg/m³)], [0], [5.3], [12], [17], [25], [63], [120], [180],
  ),
  caption: [Sorption isotherm for Cottaer Sandstone. Note the steep slope near $phi = 1$: from $phi = 0.99$ to $phi = 1.0$, water content jumps from 63 to 180~kg/m³.],
) <fig-sorption>

The derivative $d w slash d phi$ enters the capacitance matrix $bold(M)$. Near $phi = 1.0$, this derivative changes by orders of magnitude over tiny intervals of $phi$, making the system extremely stiff and difficult to solve.

== NR Algorithm

#block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%)[
  *Input:* $bold(U)_0$ = initial guess (previous timestep), $omega$ = relaxation parameter, $"tol"$ = convergence tolerance

  + Assemble $bold(A)(bold(U)_0)$, $bold(B)(bold(U)_0)$
  + For $k = 0, 1, 2, dots$
    + Compute residual: $bold(r) = bold(B) - bold(A) dot bold(U)_k$
    + Solve for correction: $bold(A) dot Delta bold(U) = bold(r)$
    + Limit increment per DOF, return `allClamped` flag (see @sec-clamping)
    + *If `allClamped` is true:* exit *without* applying $Delta bold(U)$ --- the correction is below the round-off filter threshold and applying it would seed asymmetric linear-solver noise (see @sec-roundoff-filter)
    + Apply adaptive damping (see @sec-damping)
    + Backtracking line search (see @sec-linesearch): starting from $omega_"eff"$, find the largest step that reduces the residual norm. Returns the accepted trial state $bold(U)_(k+1)$ and the reassembled $bold(A)$, $bold(B)$.
    + Check convergence (augmented metric, see @sec-roundoff-filter):
      - $"normMetric" = abs(norm(bold(U)_(k+1)) - norm(bold(U)_k)) / norm(bold(U)_(k+1))$
      - $"componentMetric" = (max_i abs(Delta U_i dot omega_"eff")) / (max_i abs(U_(k+1,i)))$
      - $"metric" = max("normMetric", "componentMetric")$
    + If $"metric" <= "tol"$: converged, exit
    + If oscillation detected: apply midpoint averaging (see @sec-oscillation), exit
    + If $k >$ MaxIterations: not converged, exit
]

The convergence criterion combines two checks: a *global* measure (relative change in the $L^2$ norm of the solution between consecutive iterations) and a *local* measure (the largest applied per-DOF correction relative to the largest solution component). Both must be below `tol` for convergence to be declared. The global metric alone has a null space --- mean-preserving (anti-symmetric) perturbations leave the $L^2$ norm unchanged --- so it can falsely report convergence after a single iteration that injects high-frequency noise. The per-component metric closes that hole. See @sec-roundoff-filter for the worked example that motivated the augmentation.

== Linear Problems

When both boundary conditions and element properties are linear (no state-dependent materials), the system $bold(A) dot bold(U) = bold(B)$ is solved in a single direct step without iteration. The solver detects this automatically by checking whether any element or boundary condition reports nonlinear behavior.

// ═════════════════════════════════════════════════════════
= Adaptive Time-Accumulation with Cross-Coupling (Outer Level)
// ═════════════════════════════════════════════════════════

The coupled hygrothermal problem is solved using a *partitioned scheme with adaptive time-stepping at the multi-domain level*. Rather than assembling a monolithic system, each domain is solved independently and then cross-domain data is exchanged after every successful sub-step.

== Design Rationale: Why Not a Staggered Loop?

The natural approach to solving coupled equations is a *staggered (Gauss--Seidel) iteration*: solve the moisture equation with frozen temperature, then solve the thermal equation with frozen humidity, repeat until both converge. Earlier versions of HygroThermFEM used this approach.

The problem emerged with extreme conditions (high humidity at elevated temperatures). Consider the `ExtremeHumidityAndTemperature` test case:

#block(fill: luma(235), inset: 12pt, radius: 4pt, width: 100%)[
  *Test setup:* Cottaer Sandstone slab, 3 columns of nodes (6 nodes total), initial $T = 80 °C$, $phi = 0.9999$, boundary condition: $T_"air" = 20 °C$, $h_c = 10$~W/(m²K), $phi_"air" = 1.0$, timestep $Delta t = 3600$~s.
]

With the old staggered approach:
- The moisture domain could not converge at $Delta t = 3600$~s (the sorption curve is extremely steep near $phi = 1.0$).
- It subdivided to $Delta t = 360$~s, then $Delta t = 36$~s, then $Delta t = 3.6$~s.
- At $Delta t = 3.6$~s, it needed $3600 slash 3.6 = 1000$ sub-steps.
- Each sub-step required multiple NR iterations.
- *All 40,000 sub-steps ran with the temperature frozen at 80~°C* --- no coupling exchange occurred during the internal grind.
- Only after completing all sub-steps did the staggered loop update the temperature, requiring the whole process to repeat.
- Total: *~15.5 seconds* for 2 timesteps, with poor coupling convergence (error ~0.01).

== The Current Approach: Caller-Level Time Accumulation

The key insight: *move the time-accumulation loop out of the individual domains and into `MultiDomain::transient`*. Each domain returns after its first successful NR solve (at whatever timestep it could converge), and cross-coupling data is exchanged after every sub-step.

== Transient Coupling Algorithm

#block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%)[
  *Input:* $bold(T)_0$ = initial temperature, $bold(H)_0$ = initial humidity, $Delta t$ = target timestep

  + $t_"total" = 0$
  + While $t_"total" < Delta t$:
    + $Delta t_"eff" = Delta t - t_"total"$ #h(2em) _(the remaining time to cover)_

    + Call moisture domain with $Delta t_"eff"$:
      - Inside the domain, NR is attempted at $Delta t_"eff"$.
      - If NR converges: domain returns the solution and $Delta t_"eff"$ unchanged.
      - If NR fails: domain divides $Delta t_"eff"$ by 10, tries again, repeats up to 3 times. Returns the solution at whatever smaller $Delta t$ succeeded. We call this $Delta t_"moisture"$.
      - Update: $Delta t_"eff" = Delta t_"moisture"$ #h(1em) _(may be smaller than what we asked for)_

    + Call thermal domain with $Delta t_"eff"$ (the value moisture actually used):
      - Same probing logic. If thermal converges at $Delta t_"eff"$: done.
      - If thermal needs an even smaller $Delta t_"thermal"$:
        - Update $Delta t_"eff" = Delta t_"thermal"$
        - Discard the moisture solution (it was for a larger $Delta t$)
        - Re-call moisture at $Delta t_"eff"$ so both domains use the same $Delta t$

    + Both domains have now solved at the same $Delta t_"eff"$. Accept results:
      - $bold(T)_"current" = bold(T)_"new"$, #h(1em) $bold(H)_"current" = bold(H)_"new"$
      - Update all nodes with both fields (cross-coupling exchange)
      - $t_"total" = t_"total" + Delta t_"eff"$

  + Coupling convergence probe (see @sec-convergence-probe)
  + Final update: advance "previous timestep" values in nodes
]

The key point: the caller asks for $Delta t_"eff"$ but the domain may return a smaller value. The caller respects whatever the domain could actually handle, ensures both domains use the same $Delta t$, and keeps looping until the full target $Delta t$ is covered.

== Walkthrough: How Sub-Steps Accumulate

To make the algorithm concrete, let us trace through the first few sub-steps for a case where the caller requests $Delta t = 3600$~s.

=== Sub-step 1 (starting from $t_"total" = 0$)

The remaining time is $3600 - 0 = 3600$~s, so we set $Delta t_"eff" = 3600$~s.

*Moisture domain is called with $Delta t = 3600$~s:*
- Inside `IDomain::transient`, the moisture domain tries NR at $Delta t = 3600$~s. The sorption curve is too steep --- NR does not converge in 50 iterations.
- The domain divides: $3600 / 10 = 360$~s. Tries NR at $Delta t = 360$~s --- converges.
- Returns the solution and reports $Delta t_"actual" = 360$~s.

Now $Delta t_"eff" = min(3600, 360) = 360$~s.

*Thermal domain is called with $Delta t = 360$~s:*
- Thermal tries NR at $Delta t = 360$~s --- converges (thermal is typically easier).
- Returns the solution and reports $Delta t_"actual" = 360$~s.

*Accept and exchange:* Both solutions are accepted. Nodes are updated with the new temperature *and* new humidity. $t_"total" = 0 + 360 = 360$~s.

=== Sub-step 2 (starting from $t_"total" = 360$~s)

Remaining time: $3600 - 360 = 3240$~s. Set $Delta t_"eff" = 3240$~s.

*Moisture domain is called with $Delta t = 3240$~s:*
- Tries $3240$~s --- fails. Divides to $324$~s --- converges.
- Returns $Delta t_"actual" = 324$~s.

$Delta t_"eff" = min(3240, 324) = 324$~s.

*Thermal domain is called with $Delta t = 324$~s:*
- Converges. Returns $Delta t_"actual" = 324$~s.

*Accept and exchange.* $t_"total" = 360 + 324 = 684$~s.

=== Sub-steps 3, 4, ... (pattern continues)

Each sub-step:
+ Asks moisture to solve for the remaining time.
+ Moisture probes and finds the largest $Delta t$ it can handle.
+ Thermal is solved at that same $Delta t$.
+ Both solutions are accepted, coupling data exchanged.
+ $t_"total"$ advances.

Note that the effective $Delta t$ can *change* from sub-step to sub-step --- it depends on what the moisture domain can handle at its current state. Early sub-steps may need smaller $Delta t$ (steep sorption region), while later sub-steps may converge at larger $Delta t$ as the solution moves to a smoother part of the curve.

=== When $t_"total"$ reaches $3600$~s

The loop exits. Typically this takes 10--20 sub-steps for moderate cases (where moisture converges at $Delta t approx 200$--$400$~s) or 50--70 sub-steps for extreme cases (where the smallest successful $Delta t$ is 3.6~s or less).

=== The key difference from the old approach

In the old solver, all sub-stepping happened *inside* the moisture domain. The thermal domain was called only once, after all moisture sub-steps completed. This meant:
- 10 moisture sub-steps at 360~s each $arrow$ all 10 used the *same* frozen temperature field
- The thermal domain saw the final moisture result but had no influence during the sub-stepping

With the new approach, thermal is solved and coupling data exchanged after *every* sub-step. Sub-step 2 sees the temperature from sub-step 1, sub-step 3 sees the temperature from sub-step 2, and so on.

== What Happens When Domains Need Different Timesteps

Sometimes the moisture domain converges at a certain $Delta t$ but the thermal domain needs a smaller one (or vice versa). The algorithm handles this through *re-solving*:

#block(fill: luma(245), inset: 10pt, radius: 3pt, width: 100%)[
  *Scenario:* Moisture converges at $Delta t = 360$~s, but thermal can only converge at $Delta t = 36$~s.

  + Moisture is called with $Delta t_"eff" = 3600$~s (remaining time).
  + Moisture probes: 3600~s fails, 360~s succeeds $arrow$ returns solution at $Delta t = 360$~s.
  + $Delta t_"eff"$ is reduced to 360~s.
  + Thermal is called with $Delta t_"eff" = 360$~s.
  + Thermal probes: 360~s fails, 36~s succeeds $arrow$ returns solution at $Delta t = 36$~s.
  + Since $Delta t_"thermal" = 36 < Delta t_"eff" = 360$: the moisture solution at 360~s is *discarded*.
  + $Delta t_"eff"$ is reduced to 36~s.
  + Moisture is *re-solved* at $Delta t = 36$~s (this time it converges easily since 36 < 360).
  + Both solutions accepted at $Delta t = 36$~s. Exchange coupling data. Advance $t_"total"$ by 36~s.
]

This ensures both domains always advance by the *same* $Delta t$, maintaining temporal consistency. In practice, the thermal domain rarely needs a smaller timestep than moisture --- the moisture equation is typically the stiffer one due to the steep sorption isotherm.

== Example: The Extreme Test Case <sec-extreme-example>

With the new approach on the same `ExtremeHumidityAndTemperature` test:

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.5pt,
    inset: 8pt,
    [*Metric*], [*Old (staggered)*], [*New (accumulation)*],
    [Moisture sub-steps (level 1)], [400], [25],
    [Moisture sub-steps (level 2)], [40,000], [25],
    [Moisture sub-steps (level 3)], [0], [22],
    [Thermal sub-steps], [0], [0],
    [Wall-clock time (2 timesteps)], [~15.5 s], [~0.3 s],
    [Coupling error (humidity)], [$1.3 times 10^(-2)$], [$1.2 times 10^(-8)$],
    [Coupling error (temperature)], [$3.8 times 10^(-1)$], [$1.6 times 10^(-7)$],
  ),
  caption: [Performance comparison on the extreme test case ($T = 80 °C$, $phi = 0.9999$). The new approach is ~50× faster with 6--7 orders of magnitude better coupling convergence.],
) <fig-comparison>

Why such a dramatic improvement? In the old approach, the moisture domain ground through thousands of sub-steps with frozen $T = 80 °C$. Since the saturation concentration $c_"sat"(80 °C)$ is very high, the moisture equation was solving with incorrect transport coefficients the entire time. In the new approach, the thermal domain is solved and coupling data exchanged after every sub-step, so the moisture domain always sees the correct (evolving) temperature field.

== Effect on Physical Results

The different time-stepping strategies produce different numerical results. Compare the first-timestep water content for three humidity levels:

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.5pt,
    inset: 8pt,
    [*Test case*], [*$phi_0 = 0.999$, $T_0 = 0 °C$*], [*$phi_0 = 0.9999$, $T_0 = 80 °C$*],
    [Old: node 1 water content], [103.2 kg/m³], [180.0 kg/m³ (saturated)],
    [New: node 1 water content], [180.0 kg/m³ (saturated)], [84.9 kg/m³],
    [Old: node 5 water content], [103.4 kg/m³], [180.0 kg/m³ (saturated)],
    [New: node 5 water content], [27.8 kg/m³], [53.2 kg/m³],
  ),
  caption: [Water content after first timestep. The new solver produces different --- and physically more meaningful --- moisture distributions because temperature evolves concurrently with moisture.],
)

The new results reflect the fact that at $T = 80 °C$, the saturation concentration is very high, so the initial $phi = 0.9999$ does not correspond to as much liquid water as it does at $T = 0 °C$. The old solver did not capture this because the temperature was frozen during moisture sub-stepping.

== Coupling Convergence Probe <sec-convergence-probe>

After the time-accumulation loop completes, the solver performs a *coupling convergence probe*: each domain is re-solved at a tiny timestep ($Delta t_"probe" = Delta t times 10^(-6)$) from the current state, and the resulting solution change is measured. If the sub-stepping with frequent coupling did its job, the probe should show negligible change.

The probe errors are reported as `temperatureError` and `humidityError` in the `Solution` struct. Typical values:

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.5pt,
    inset: 8pt,
    [*Test case*], [*Humidity coupling error*], [*Temperature coupling error*],
    [$phi = 0.999$, $T = 0 °C$], [$1.8 times 10^(-6)$], [$4.0 times 10^(-7)$],
    [$phi = 0.9999$, $T = 30 °C$], [$2.7 times 10^(-7)$], [$6.3 times 10^(-9)$],
    [$phi = 0.9999$, $T = 80 °C$], [$1.2 times 10^(-8)$], [$1.6 times 10^(-7)$],
  ),
  caption: [Coupling convergence probe errors. All are well below the NR convergence tolerance of $10^(-5)$, confirming self-consistent solutions.],
)

== Steady-State Coupling

The steady-state coupling follows a simpler pattern: alternately solve each domain and exchange data, repeating until both converge. No time-accumulation is needed since there is no time derivative.

// ═════════════════════════════════════════════════════════
= Convergence Stabilization Techniques
// ═════════════════════════════════════════════════════════

Near-saturation humidity conditions ($phi$ close to 1.0) present particular challenges for convergence. The steep sorption isotherm (@fig-sorption) means that small changes in $phi$ produce large changes in water content and transport coefficients, leading to several pathological behaviors in the NR iteration.

== Per-DOF Increment Limiting <sec-clamping>

*Problem.* The NR correction $Delta bold(U)$ can project the humidity solution outside the physical range $[0, 1]$. For example, if $phi_i = 0.98$ and $Delta phi_i = 0.05$, the projected value would be $phi_i = 1.03$ --- physically meaningless. Post-process clamping (forcing $phi = 1.0$) introduces an inconsistency: the system matrices were assembled assuming the unconstrained state, but the actual solution is clamped. On the next iteration, the matrices reflect $phi = 1.0$ (with its enormous $d w slash d phi$), but the correction was computed for $phi = 1.03$. This causes the solver to repeatedly overshoot and clamp, creating artificial oscillations.

*Solution.* Before applying the correction, each DOF's increment is scaled so the projected value $U_i + omega dot Delta U_i$ remains within $[0, 1]$:

$ "For each DOF" i: quad cases(
  Delta U_i = (1 - U_i) / omega & "if" U_i + Delta U_i dot omega > 1,
  Delta U_i = -U_i / omega & "if" U_i + Delta U_i dot omega < 0,
) $

This keeps the solution within physical bounds at every iteration while maintaining consistency between the solution state and the assembled system matrices.

*Early exit on negligible correction (round-off filter).* `limitIncrement` also returns a flag indicating whether *every* DOF correction has $|Delta U_i dot omega| < 10^(-12)$ after clamping. When true, two situations are possible: either the solution genuinely sits at a fixed point of the iteration (e.g., pinned at $phi = 0$ or $phi = 1$ at every node, or at a steady state interior to the bounds), *or* the residual $bold(r) = bold(B) - bold(A) bold(U)$ is at the floating-point noise floor and the linear solver is returning round-off as $Delta bold(U)$. Either way, applying the correction is wrong --- in the first case it changes nothing meaningful, in the second case it injects asymmetric noise that the norm-based convergence check cannot detect (see @sec-roundoff-filter). The NR loop therefore exits *without applying* $Delta bold(U)$ when this flag fires, treating the current iterate as the solution.

=== Example

In the `ExtremeHumidityAndTemperature` test, after the first NR iteration of the moisture domain, several nodes overshoot $phi = 1.0$. Without increment limiting, the solver would oscillate between $phi = 0.97$ and $phi = 1.03$ indefinitely. With limiting, the overshooting DOFs are gently clamped to $phi = 1.0$ on each iteration until the system naturally settles at saturation.

== Adaptive Damping <sec-damping>

*Problem.* In strongly nonlinear problems, the NR correction norm $norm(Delta bold(U))$ can grow explosively between iterations rather than decreasing. This indicates the solver is diverging --- each correction overshoots further than the previous one.

*Solution.* Track the correction norm across iterations. When it grows by a factor exceeding 1000 on two or more consecutive iterations (after the solver has had at least $"MaxIterations" slash 2$ warmup iterations), the damping factor is halved:

$ omega_"eff" = omega dot alpha, quad "where" alpha in [1 slash 16, 1] $

The damping factor $alpha$ starts at 1.0 (no damping). When explosive growth is detected, $alpha$ is halved (minimum $1 slash 16$). When growth stops, $alpha$ is doubled back toward 1.0. This ensures that normal NR fluctuations are not penalized, while genuinely divergent behavior is controlled.

The warmup period ($"MaxIterations" slash 2$ iterations before damping activates) prevents interference with the early NR phase, where corrections naturally fluctuate in strongly nonlinear regions before the iteration finds its basin of convergence.

== Oscillation Detection and Midpoint Averaging <sec-oscillation>

*Problem.* Piecewise-linear material data (e.g., the sorption isotherm in @fig-sorption) can have steep kinks where the derivative $d w slash d phi$ changes abruptly. When the NR iteration crosses such a kink, the Jacobian (effectively the assembled $bold(A)$ matrix) changes discontinuously, causing the solution to alternate between two states on consecutive iterations (a "2-cycle"). The convergence metric stabilizes at a constant nonzero value instead of decreasing toward zero.

*Detection.* After a minimum warmup period (4 iterations), the solver compares the convergence metric between consecutive iterations. If the relative change is less than 1%:

$ abs("metric"_k - "metric"_(k-1)) / "metric"_(k-1) < 0.01 $

the solver concludes that a 2-cycle is present.

*Resolution.* The solution is set to the midpoint of the last two iterates:

$ bold(U)_"final" = 1/2 (bold(U)_k + bold(U)_(k-1)) $

This places the solution at the center of the oscillation band, which lies at (or very near) the material property kink. The midpoint is the best physically meaningful estimate when the true solution sits exactly on a non-smooth point of the material curve.

=== Why 2-Cycles Occur at Sorption Kinks

Consider a node where $phi$ oscillates between 0.989 and 0.991 across iterations:
- At $phi = 0.989$: the sorption derivative $d w slash d phi approx 600$ (below the kink at $phi = 0.99$, where $w$ jumps from 17 to 63 kg/m³ over $Delta phi = 0.04$).
- At $phi = 0.991$: the sorption derivative $d w slash d phi approx 4000$ (above the kink).
- The capacitance matrix $bold(M)$ changes by nearly an order of magnitude between the two states, so the system dynamics are fundamentally different.
- The solution computed with the "low" Jacobian overshoots into the "high" region, and vice versa.

The midpoint $phi = 0.990$ places the solution right at the kink, which is consistent with both the low and high material data.

== Backtracking Line Search <sec-linesearch>

*Problem.* Newton--Raphson can overshoot wildly in stiff regions where the linearization badly misrepresents the true nonlinear behavior. This is most pronounced for the moisture equation near sorption-curve saturation, where $d w slash d phi$ can change by orders of magnitude across a small $phi$ window. Without line search, the iteration can enter a deterministic multi-cycle that bounces between a near-uniform state and a saturated state, accept a non-physical midpoint via the oscillation handler, and commit a mass-non-conserving result.

*Solution.* After computing the damped correction $omega_"eff" dot Delta bold(U)$, the solver performs an Armijo-style backtracking line search. The full damped step is tried first; if it makes the residual grow, the step length is halved and the system is reassembled at the smaller trial state. This repeats for up to 8 attempts:

$ "For" ell = 0, 1, dots, 7: $
$ quad bold(U)_"trial" = bold(U)_k + omega_"eff" dot 2^(-ell) dot Delta bold(U) $
$ quad "Reassemble" bold(A)(bold(U)_"trial")", " bold(B)(bold(U)_"trial")) $
$ quad "If" norm(bold(B) - bold(A) bold(U)_"trial") < norm(bold(r)_k): "accept" bold(U)_"trial"", exit" $

If no attempt reduces the residual, the last attempt ($ell = 7$, step reduced by $1 slash 128$) is accepted regardless. This ensures the iteration always makes progress, even if the full step is severely overshooting.

The line search adds cost (each backtrack requires a reassembly), but it only fires when the residual would otherwise grow. In well-behaved regions, the full step is accepted on the first attempt with zero overhead beyond one residual norm evaluation.

*Implementation note.* The line search is encapsulated in `IDomain::backtrackingLineSearch`, which returns a `LineSearchResult` struct containing the accepted trial solution, the reassembled system matrices, and the effective relaxation factor used. This keeps the NR loop body (`performNRIteration`) focused on orchestration.

== Smooth Vapor Transfer Coefficient

*Problem.* The water vapor transfer coefficient $beta$ originally had a hard step discontinuity at $phi = 1.0$:

$ beta = cases(beta_"value" & "if" phi <= 1.0, 0 & "if" phi > 1.0) $

During NR iteration, small overshoots past $phi = 1.0$ would zero out the transfer coefficient, then on the next iteration (after clamping back to 1.0) it would jump back to full value. This created a Jacobian discontinuity at the exact point where the solver was trying to converge.

*Solution.* Replace the step function with a smooth linear ramp:

$ "smoothFactor" = "clamp"((1 + w - phi) / w, space 0, space 1), quad w = 0.01 $
$ beta = beta_"value" dot "smoothFactor" $

This produces: full transfer below $phi = 1.0$, a linear taper over $[1.0, 1.01]$, and zero transfer above $phi = 1.01$. The narrow transition width ensures negligible impact on physical results while eliminating the Jacobian discontinuity.

== Round-Off Filter and Augmented Convergence Metric <sec-roundoff-filter>

*Problem.* Two interacting defects in the original NR loop allowed pure floating-point round-off to grow into a catastrophic numerical drift on near-uniform fields. The defects only showed up on problems where:

+ The sorption derivative $d w slash d phi$ is large (high humidity near $phi = 1$).
+ The mass term $bold(M)$ dominates the conductance $bold(K)$ in $bold(A) = bold(M) slash Delta t + bold(K) + bold(H)$ at the requested $Delta t$.
+ The expected solution is essentially uniform in space (matching boundary conditions, no thermal gradient), so the true residual is zero up to round-off.

#block(fill: luma(235), inset: 12pt, radius: 4pt, width: 100%)[
  *Worked example: stucco beam at $phi = 0.999$, $T = 30 °C$, matching boundary conditions on both edges, $Delta t = 3600$~s, moisture domain only.* The expected behavior is a perfectly stationary solution. The original solver instead produced a slow drift over the first six timesteps that suddenly collapsed at step 7, with the level-1 NR failing to converge and the subdivider falling back to $Delta t = 360$~s --- which converged to a solution with water content dropped from 110~kg/m³ to ~93~kg/m³ in a single sub-step (a 17% jump on a problem where the answer should never move).
]

*Defect 1 --- Round-off applied as a real correction.* At iteration 0 of timestep 0, the residual $bold(r) = bold(B) - bold(A) bold(U)_0$ was at the floating-point noise floor ($~10^(-19)$). The linear solver, given a near-zero RHS, returned $Delta bold(U) tilde.op 10^(-16)$ --- the LU factorization's intrinsic round-off. *Crucially, the round-off was asymmetric:* the magnitudes at the eight nodes of the test mesh differed by an order of magnitude, with no left/right symmetry. The original NR loop applied this $Delta bold(U)$ unconditionally, *then* checked the `allClamped` flag (which was true, since every component was below $10^(-12)$) and exited. The result: the perfectly uniform initial state $bold(U)_0$ became a noisy state $bold(U)_1 = bold(U)_0 + Delta bold(U)_"noise"$, with the noise pattern set by the linear solver's internal LU rounding decisions, not by physics.

*Defect 2 --- Norm-based convergence metric is blind to the noise.* The original convergence criterion was

$ "metric" = abs(norm(bold(U)_(k+1)) - norm(bold(U)_k)) / norm(bold(U)_(k+1)) <= "tol" $

This $L^2$ norm difference has a *null space*: any perturbation $delta bold(U)$ that is orthogonal to $bold(U)$ (or, more loosely, any zero-mean perturbation) leaves the norm essentially unchanged. The asymmetric noise injected by Defect 1 was exactly such a perturbation. So even after the noise had grown to $tilde.op 10^(-7)$ over six timesteps and started producing a *real* residual ($tilde.op 10^(-8)$), and the linear solver started returning a non-trivial $Delta bold(U) tilde.op 10^(-5)$, NR would still exit after one iteration because the $L^2$ norm of the solution barely moved.

But $Delta bold(U) tilde.op 10^(-5)$ is no longer round-off. It is what the diagonally-dominant explicit-Euler-like step computes when $bold(M) slash Delta t$ dominates $bold(A)$:

$ Delta bold(U) approx (bold(M) slash Delta t)^(-1) bold(r) approx Delta t dot bold(M)^(-1) (-bold(K) bold(U)_"noisy") $

Explicit Euler on a parabolic stencil at large $Delta t$ is unstable for high-frequency modes. Each timestep injected one explicit step's worth of noise growth and called itself converged. After seven timesteps the noise was large enough that NR overshot into physical bounds, the metric finally fired "non-converged", the subdivider kicked in, and the catastrophe became visible in the water content.

*Why the steep sorption curve matters.* Both defects are amplified by $d w slash d phi$ being huge near $phi = 1$:

+ It inflates the lumped mass $bold(M) = rho dot (d w slash d phi) dot V$, which makes $bold(M)$ dominate $bold(A)$, which turns NR's first step into the unstable explicit-Euler step described above.
+ It is also the conversion factor from $phi$-noise to water-content-noise, so a $10^(-16)$ perturbation in $phi$ becomes a visible $tilde.op 10^(-13)$ perturbation in $w$, then grows from there.

A flat sorption curve hides the bug entirely. The high-humidity regime exposes both flaws simultaneously.

*Solution part 1 --- Round-off filter (Defect 1).* When `limitIncrement` reports `allClamped = true`, exit *before* applying $Delta bold(U)$. See the early-exit step in the @sec-clamping NR pseudocode. The current iterate is treated as the converged result, leaving $bold(U)$ untouched. This stops the noise from being seeded in the first place.

*Solution part 2 --- Per-component convergence metric (Defect 2).* Augment the convergence check with a per-component measure that has no null space:

$ "componentMetric" = (max_i abs(Delta U_i dot omega_"eff")) / (max_i abs(U_(k+1,i))) $
$ "metric" = max("normMetric", "componentMetric") $

The combined metric is the strictly more conservative of the two. Convergence requires that *both* the global norm has stagnated *and* every per-DOF move is small relative to the largest solution component. For the worked example above, the per-component metric correctly registers the $10^(-5)$ move at step 6 as non-trivial relative to $|U_i| approx 1$, forcing NR to keep iterating instead of exiting on the loose norm criterion.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, center, center),
    stroke: 0.5pt,
    inset: 8pt,
    [*Quantity*], [*Before fix*], [*After fix*],
    [Water content over 10 steps], [drift then collapse to ~89 kg/m³], [exactly 110.000000 kg/m³],
    [`waterContentError` at step 9], [$2.5 times 10^(-9)$], [$0$],
    [`progressMoisture` (level-1, level-2, level-3)], [(17, 2, 0)], [(0, 0, 0)],
    [Total NR calls in moisture trace], [~580 (catastrophic cascade at step 7)], [20 (one per main call + one per probe)],
  ),
  caption: [Effect of the round-off filter and augmented convergence metric on the stucco-beam uniform-field test.],
)

*Interaction with the existing fixes.* The round-off filter and the per-component metric do *not* override the increment-limit clamping (@sec-clamping), the adaptive damping (@sec-damping), or the oscillation midpoint (@sec-oscillation). Those mechanisms still fire in their respective regimes; the new checks add a noise-floor sanity gate at each end of the NR iteration --- "do not move if the move is round-off" before applying, and "do not declare converged if any individual component is still moving" after applying.

// ═════════════════════════════════════════════════════════
= Adaptive Timestep Probing (Middle Level)
// ═════════════════════════════════════════════════════════

When the NR iteration fails to converge at a requested timestep, the domain does not give up. Instead, it *subdivides* and tries progressively smaller timesteps until one succeeds, then *returns immediately* with the first successful solution and the actual $Delta t$ used.

== Algorithm

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

== Probing Depths

With default settings (3 levels, 10 subdivisions per level), the probing sequence for a requested $Delta t = 3600$~s is:

#figure(
  table(
    columns: (auto, auto, auto),
    align: (center, center, left),
    stroke: 0.5pt,
    inset: 8pt,
    [*Level*], [*$Delta t_"probe"$*], [*Comment*],
    [0], [3600 s], [Try full timestep first],
    [1], [360 s], [First subdivision],
    [2], [36 s], [Second subdivision],
    [3], [3.6 s], [Third subdivision (last attempt)],
  ),
  caption: [Timestep probing sequence. Each level divides by 10. If level 3 fails, the solver throws an error.],
)

The number of sub-steps needed to cover the target $Delta t$ depends on which level succeeds. If level 2 succeeds ($Delta t = 36$~s), the outer loop needs $3600 slash 36 = 100$ sub-steps. Each sub-step exchanges cross-coupling data, so the moisture field sees 100 updates to the temperature field over the course of one target timestep.

== Observer Notifications

Each probing level is reported to registered `TimestepObserver` instances. Tests use this to verify the expected probing depth:

```cpp
TestHelper::ObserveSimulationProgress progressMoisture;
multiDomain.subscribeMoisture(&progressMoisture);
// ... run simulation ...
// calls() returns a vector where index 0 = level-1 count, etc.
EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));  // no subdivisions needed
```

// ═════════════════════════════════════════════════════════
= Worked Example: Three Humidity Regimes
// ═════════════════════════════════════════════════════════

The `MultiDomain_HighHumidity` test suite exercises the solver at three progressively harder humidity levels, all using the same Cottaer Sandstone material and a simple 2-element mesh (6 nodes).

== Mesh Geometry

The test mesh consists of 6 nodes arranged in a 3$times$2 grid, forming 2 quadrilateral elements. The x-coordinates are $x in {0.15, 0.05, 0.00}$~m and y-coordinates are $y in {0.00, 0.05}$~m. The convective boundary condition is applied on the right edge (nodes 5--6).

#figure(
  cetz.canvas(length: 1cm, {
    import cetz.draw: *

    // Scale: 1cm = 0.02m, so 0.15m = 7.5cm, 0.05m = 2.5cm
    let sx = 50   // scale x: multiply meters by this
    let sy = 50   // scale y

    let xs = (0.15, 0.05, 0.00)
    let ys = (0.00, 0.05)

    // Draw elements (filled)
    rect((xs.at(1)*sx, ys.at(0)*sy), (xs.at(0)*sx, ys.at(1)*sy),
         fill: luma(240), stroke: 0.5pt, name: "e1")
    rect((xs.at(2)*sx, ys.at(0)*sy), (xs.at(1)*sx, ys.at(1)*sy),
         fill: luma(230), stroke: 0.5pt, name: "e2")

    // Element labels
    content((3.75, 1.25), text(size: 9pt, fill: luma(120))[Element 1])
    content((1.25, 1.25), text(size: 9pt, fill: luma(120))[Element 2])

    // Node dots and labels
    let nodes = (
      (xs.at(0), ys.at(0), "1"),
      (xs.at(0), ys.at(1), "2"),
      (xs.at(1), ys.at(0), "3"),
      (xs.at(1), ys.at(1), "4"),
      (xs.at(2), ys.at(0), "5"),
      (xs.at(2), ys.at(1), "6"),
    )

    for (nx, ny, label) in nodes {
      circle((nx*sx, ny*sy), radius: 0.12, fill: black)
      content((nx*sx, ny*sy - 0.45),
              text(size: 8pt, weight: "bold")[Node #label])
    }

    // Boundary condition arrows on right edge
    for yoff in (0.5, 1.0, 1.5, 2.0) {
      line((-0.8, yoff), (-0.15, yoff),
           stroke: (paint: rgb("#c44"), thickness: 1.2pt),
           mark: (end: ">", fill: rgb("#c44")))
    }
    content((-1.8, 1.25),
            text(size: 8pt, fill: rgb("#c44"))[BC: $T_"air"$, $phi_"air"$\ $h_c = 10$])

    // Dimension annotations
    set-style(stroke: (paint: luma(150), thickness: 0.5pt, dash: "dashed"))
    line((0, -0.6), (2.5, -0.6), mark: (start: "|", end: "|"))
    content((1.25, -1.0), text(size: 7pt, fill: luma(100))[0.10 m])
    line((2.5, -0.6), (7.5, -0.6), mark: (start: "|", end: "|"))
    content((5.0, -1.0), text(size: 7pt, fill: luma(100))[0.10 m])
    line((7.8, 0), (7.8, 2.5), mark: (start: "|", end: "|"))
    content((8.6, 1.25), text(size: 7pt, fill: luma(100))[0.05 m])
  }),
  caption: [Two-element mesh used in the `MultiDomain_HighHumidity` tests. Red arrows indicate the convective boundary condition applied at nodes 5 and 6. All nodes share the same initial state ($T_0$, $phi_0$).],
) <fig-mesh>

== Test 1: $phi = 0.999$, $T = 0 °C$ (Hard)

#block(fill: luma(235), inset: 12pt, radius: 4pt, width: 100%)[
  *Setup:* All nodes start at $T = 0 °C$, $phi = 0.999$. Boundary condition at nodes 5--6: $T_"air" = 20 °C$, $h_c = 10$, $phi_"air" = 1.0$. Two timesteps of $Delta t = 3600$~s.
]

At $phi = 0.999$, the sorption curve gives $w = 120$~kg/m³. The boundary drives humidity toward 1.0, where $w = 180$~kg/m³. This means the moisture equation must push through the steep part of the sorption curve.

*Results after timestep 1:*
- Nodes 1--4 saturate ($w = 180$~kg/m³) --- the boundary humidity of 1.0 drives them to full saturation.
- Nodes 5--6 remain at $w approx 28$~kg/m³ --- they are farther from the boundary and lag behind.
- Temperature drops slightly below 0°C at boundary nodes (evaporative cooling effect).

*Solver effort:* 8 level-1 probes, 7 level-2 probes, 3 level-3 probes. The moisture domain needed sub-stepping, but the thermal domain converged at the full $Delta t$.

== Test 2: $phi = 0.9999$, $T = 30 °C$ (Harder)

Same mesh and boundary conditions, but starting at $T = 30 °C$ and $phi = 0.9999$. Higher temperature means higher saturation concentration --- more moisture capacity, stronger transport.

*Results after timestep 1:*
- Nodes 1--2 at $w = 81$~kg/m³ (not yet saturated).
- Nodes 3--6 saturate ($w = 180$~kg/m³).
- Temperature barely changes (30°C → ~31°C at boundary) because the thermal boundary is also at 20°C, creating a cooling gradient.

*Solver effort:* 4 level-1, 4 level-2, 1 level-3 probes. Easier than Test 1 because the higher temperature provides more "room" in the sorption curve.

== Test 3: $phi = 0.9999$, $T = 80 °C$ (Extreme) <sec-test3>

The stress test. At 80°C, the saturation concentration $c_"sat" approx 290$~g/m³ (vs. ~5~g/m³ at 0°C). The moisture transport coefficients are dramatically different.

*Results after timestep 1:*
- Water content ranges from 53 to 85~kg/m³ --- no nodes saturate because at 80°C the high saturation concentration means $phi = 0.9999$ corresponds to much less liquid water.
- Temperature drops sharply from 80°C to 6--72°C across the slab (strong evaporative cooling from the rapid moisture transport).

*Results after timestep 2:*
- All nodes saturate ($w = 180$~kg/m³) and temperature equilibrates toward 33--66°C.

*Solver effort:* 25 level-1, 25 level-2, 22 level-3 probes. This was the case that took 40,000 level-2 probes with the old staggered approach.

// ═════════════════════════════════════════════════════════
= Diagnostic Output
// ═════════════════════════════════════════════════════════

Both `IDomain` and `MultiDomain` support an optional diagnostic stream via `setDiagnosticStream(std::ostream*)`. When non-null, the solver writes CSV-formatted data for every NR iteration. When the pointer is null (the default), there is zero overhead --- no checks, no formatting, no output.

== CSV Format

*Data rows* contain one record per NR iteration:

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    stroke: 0.5pt,
    inset: 6pt,
    [`subtimestep_dt`], [The sub-timestep size being attempted],
    [`nr_iter`], [NR iteration number within this sub-step],
    [`metric`], [Convergence metric (relative norm change)],
    [`converge_tol`], [Convergence tolerance ($10^(-5)$ by default)],
    [`correction_norm`], [Raw $norm(Delta bold(U))$ before clamping/damping],
    [`damping_factor`], [Current adaptive damping factor $alpha in [1 slash 16, 1]$],
    [`all_clamped`], [1 if all DOFs hit physical bounds, 0 otherwise],
    [`converged`], [1 if this iteration achieved convergence],
    [`reason`], [Why converged: `metric`, `clamped`, `oscillation`, or empty],
    [`u0, u1, ...`], [Solution values at each DOF],
  ),
  caption: [CSV column definitions],
)

*Comment lines* (`#`) mark structural boundaries:

```
# substep=0 totalTime=0.000000e+00 remaining=3.600000e+03
# MOISTURE substep=0
subtimestep_dt,nr_iter,metric,...
3.60000000e+03,1,2.47240431e-03,...
3.60000000e+03,2,2.11268134e-02,...
# THERMAL substep=0 dt=3.600000e+03
...
# accepted dt=3.600000e+03 totalTime=3.600000e+03 ...
# CONVERGENCE CHECK
# MOISTURE convergence probe
# THERMAL convergence probe
# FINAL temp_err=3.990158e-07 hum_err=1.754033e-06
```

== Usage

```cpp
#include <fstream>
// Enable diagnostics
std::ofstream diagFile("solver_diagnostics.csv");
multiDomain.setDiagnosticStream(&diagFile);
// ... run simulation ...
// Disable diagnostics
multiDomain.setDiagnosticStream(nullptr);
```

Load in Python for analysis:
```python
import pandas as pd
df = pd.read_csv("solver_diagnostics.csv", comment="#")
# Plot convergence metric vs. iteration
df.plot(x="nr_iter", y="metric", logy=True)
```

// ═════════════════════════════════════════════════════════
= Summary of Parameters
// ═════════════════════════════════════════════════════════

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
    [Adaptive damping growth threshold], [$1000$], [Correction norm growth factor that triggers damping],
    [Minimum damping factor], [$1 slash 16$], [Floor for the adaptive damping multiplier],
    [Oscillation check threshold], [$0.01$], [Relative metric change below which 2-cycle is detected],
    [Min iterations for oscillation check], [$4$], [Warmup before oscillation detection activates],
    [Clamp tolerance], [$10^(-12)$], [Absolute increment below which a DOF is considered clamped],
    [Line search max attempts], [$8$], [Maximum backtracking halvings before accepting the step],
    [Line search shrink factor], [$0.5$], [Step length reduction per backtracking attempt],
    [Vapor transfer transition width], [$0.01$], [RH range over which $beta$ tapers to zero],
    [Convergence probe $Delta t$ factor], [$10^(-6)$], [Fraction of target $Delta t$ used for coupling probe],
  ),
  caption: [Solver configuration parameters],
)

// ═════════════════════════════════════════════════════════
= Licensing and Prior Art
// ═════════════════════════════════════════════════════════

The solver is built from a combination of three well-known, public-domain techniques:

+ *Adaptive time-stepping* (try a $Delta t$, if NR fails subdivide and retry) --- this is the standard approach in ODE/PDE solvers since the 1960s. Examples include Runge--Kutta--Fehlberg (1969), LSODE (Hindmarsh, 1983), and SUNDIALS (Hindmarsh et al., 2005). No patents or licenses apply.

+ *Partitioned (staggered) coupling* (solve domain A, exchange data, solve domain B, repeat) --- standard multi-physics approach since Felippa & Park (1980). Used in every major commercial FEM package (ANSYS, ABAQUS, COMSOL). No patents or licenses apply.

+ *Early return from failed solve* (domain returns after its first successful sub-step instead of accumulating internally) --- a software design choice about where to place a loop. No algorithmic novelty; this is ordinary control flow.

The specific innovation in HygroThermFEM is the *combination* of these three: moving the time-accumulation loop to the caller level so that cross-coupling data is exchanged after every adaptive sub-step. This is an engineering design decision --- the same operations in a different order --- not a new algorithm. The closest published concept is "subcycling" in multi-physics simulations (Belytschko et al., 1979), where different domains advance at different rates, which has been in the open literature for over 40 years.

Individual stabilization techniques are equally standard:

- *Newton--Raphson iteration* --- Isaac Newton (1685), Joseph Raphson (1690).
- *Under-relaxation* --- classical technique, described in every numerical methods textbook.
- *Per-DOF increment limiting* --- standard in constrained nonlinear solvers.
- *Midpoint averaging for oscillation resolution* --- basic numerical stabilization.
- *Adaptive damping* --- standard divergence control in NR solvers.
- *Backtracking line search (Armijo condition)* --- Armijo (1966), described in every numerical optimization textbook.
- *Smooth transition functions* --- standard regularization technique to avoid Jacobian discontinuities.

The specific combination, tuning parameters, and implementation details are original engineering work by the HygroThermFEM authors.
