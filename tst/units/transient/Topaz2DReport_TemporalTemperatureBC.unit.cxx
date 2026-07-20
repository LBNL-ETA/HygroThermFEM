#include <array>
#include <cmath>
#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.5 -- Plane Transient Linear, Temporal Temperature Boundary Condition.
///
/// Unit slab (alpha = 1) at zero initial temperature, adiabatic at x = 0, surface
/// at x = 1 held at T0 e^(beta t) with T0 = 1, beta = 0.1. Analytical solution:
/// Carslaw & Jaeger p. 132. Expected values are the exact series at the report's
/// checkpoints (hygrothermfem_python, analytic.slab_temperature_exponential); the
/// tolerance is the measured backward-Euler discretization band at dt = 0.01 on
/// ten elements (max 8.8e-3).
///
/// The per-timestep BC vector convention: element i (0-based) is the boundary
/// value DURING step i + 1, i.e. the end-of-step value backward Euler applies at
/// t = (i + 1) dt. Per-step BC vectors advance only through the explicit
/// timestepIndex argument of MultiDomain::transient, so the steps are driven
/// manually here rather than through transientMultiStep.
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_TemporalTemperatureBC, ExponentialSurface)
{
    SCOPED_TRACE("Begin Test: unit slab, exponentially rising surface temperature.");

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
    constexpr auto growthRate = 0.1;

    std::vector<double> surfaceTemperatures;
    surfaceTemperatures.reserve(nSteps);
    for(std::size_t step = 1; step <= nSteps; ++step)
    {
        surfaceTemperatures.push_back(std::exp(growthRate * static_cast<double>(step) * dTime));
    }
    multiDomain.thermal().createBC_FixedTemperature(21, 22, surfaceTemperatures);

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

    TestHelper::CsvDump dump("topaz_temporal_temperature.csv", 11);
    for(std::size_t step = 0; step < history.size(); ++step)
    {
        dump.addRow(step + 1, TestHelper::bottomRow(history[step], 11, 2));
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.050807508, 0.265508026, 1.010050167}},
      {50, {0.644476626, 0.761372104, 1.051271096}},
      {100, {0.948353477, 0.991926384, 1.105170918}},
      {200, {1.153977435, 1.171120019, 1.221402758}},
    };

    for(const auto & [step, expected] : checkpoints)
    {
        for(std::size_t pos = 0; pos < expected.size(); ++pos)
        {
            EXPECT_NEAR(expected[pos], history[step - 1][pos * 10], 0.01)
              << "step " << step << ", position " << pos;
        }
    }
}
