#include <array>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.1 -- Plane Transient Linear, Constant Temperature Boundary Condition.
///
/// Unit slab (alpha = 1) at uniform initial temperature 1.0, adiabatic at x = 0,
/// surface at x = 1 stepped to 0. Analytical solution: odd-cosine series,
/// Carslaw & Jaeger p. 97. Expected values are the exact series at the report's
/// checkpoints (hygrothermfem_python, analytic.slab_temperature_step); the
/// tolerance is the measured backward-Euler discretization band at dt = 0.01 on
/// ten elements (max 8.8e-3).
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_TemperatureBC, TransientSeries)
{
    SCOPED_TRACE("Begin Test: unit slab, surface temperature step.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 1.0,
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

    constexpr auto tSurface = 0.0;
    multiDomain.thermal().createBC_FixedTemperature(21, 22, tSurface);

    constexpr auto dTime = 0.01;
    constexpr auto nSteps = 200;

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

    for(std::size_t step = 0; step < solution.size(); ++step)
    {
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.949305363, 0.735651315, 0.0}},
      {50, {0.370777430, 0.262188276, 0.0}},
      {100, {0.107977044, 0.076351300, 0.0}},
      {200, {0.009156990, 0.006474970, 0.0}},
    };

    for(const auto & [step, expected] : checkpoints)
    {
        for(std::size_t pos = 0; pos < expected.size(); ++pos)
        {
            EXPECT_NEAR(expected[pos], solution[step - 1][pos * 10], 0.01)
              << "step " << step << ", position " << pos;
        }
    }
}
