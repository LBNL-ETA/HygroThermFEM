#include <array>
#include <cmath>
#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.6 -- Plane Transient Linear, Temporal Flux Boundary Condition.
///
/// Unit slab (alpha = 1, k = 1) at zero initial temperature, adiabatic at x = 0,
/// surface heat flux q = q0 / sqrt(t) at x = 1 (the m = -1 member of the report's
/// q0 t^(m/2) family). Analytical solution: image-source erfc series, Carslaw &
/// Jaeger p. 113. Expected values are the exact series at the report's
/// checkpoints (hygrothermfem_python, analytic.slab_flux_inverse_sqrt); the
/// tolerance is the measured backward-Euler discretization band at dt = 0.01 on
/// ten elements (max 1.7e-2, largest at the singular start).
///
/// The point flux is singular at t = 0, so each per-timestep value is the
/// STEP-AVERAGE of q over [t_(i-1), t_i], q_i = 2 q0 (sqrt(t_i) - sqrt(t_(i-1)))
/// / dt, which preserves the injected-heat integral. End-of-step point sampling
/// would leave a persistent deficit of about 1.46 q0 sqrt(dt). The reference
/// solver's power_average flux profile applies the same rule. Per-step BC
/// vectors advance only through the explicit timestepIndex argument of
/// MultiDomain::transient, so the steps are driven manually.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_TemporalFluxBC, InverseSqrtFlux)
{
    SCOPED_TRACE("Begin Test: unit slab, q0 / sqrt(t) surface flux.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

    auto params = TestHelper::TestMaterial();
    params.density = 1.0;
    params.heatCapacity = 1.0;
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0})
        .height(0.05)
        .material(material.name())
        .state(state)
        .build();

    constexpr auto dTime = 0.01;
    constexpr std::size_t nSteps = 200;
    constexpr auto fluxCoefficient = 1.0;

    // Step-average of q0 / sqrt(t); positive heats the slab on this segment.
    std::vector<double> surfaceFluxes;
    surfaceFluxes.reserve(nSteps);
    for(std::size_t step = 1; step <= nSteps; ++step)
    {
        const double tNew = static_cast<double>(step) * dTime;
        const double tOld = static_cast<double>(step - 1) * dTime;
        surfaceFluxes.push_back(
          2.0 * fluxCoefficient * (std::sqrt(tNew) - std::sqrt(tOld)) / dTime);
    }
    multiDomain.thermal().createBC_FixedFlux(21, 22, surfaceFluxes);

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    auto humidities = multiDomain.nodes().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> history;
    history.reserve(nSteps);
    for(std::size_t step = 0; step < nSteps; ++step)
    {
        const auto stepResult = multiDomain.transient(temperatures, humidities, dTime, step);
        temperatures = stepResult.temperature;
        humidities = stepResult.humidity;
        history.push_back(temperatures);
    }

    TestHelper::CsvDump dump("topaz_temporal_flux.csv", 11);
    for(std::size_t step = 0; step < history.size(); ++step)
    {
        dump.addRow(step + 1, TestHelper::bottomRow(history[step], 11, 2));
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.089853905, 0.468545925, 1.772481303}},
      {50, {1.134409023, 1.353411911, 1.933972637}},
      {100, {1.821383046, 1.957699975, 2.346725524}},
      {200, {2.706823821, 2.798745120, 3.068381660}},
    };

    for(const auto & [step, expected] : checkpoints)
    {
        for(std::size_t pos = 0; pos < expected.size(); ++pos)
        {
            EXPECT_NEAR(expected[pos], history[step - 1][pos * 10], 0.02)
              << "step " << step << ", position " << pos;
        }
    }
}
