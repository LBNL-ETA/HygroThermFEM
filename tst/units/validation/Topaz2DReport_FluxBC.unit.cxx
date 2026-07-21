#include <array>
#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.2 -- Plane Transient Linear, Constant Flux Boundary Condition.
///
/// Unit slab (alpha = 1, k = 1) at uniform initial temperature 0, adiabatic at
/// x = 0, uniform surface heat flux q = 1 into the slab at x = 1. The engine's
/// flux sign at a boundary segment follows the segment's node ordering: on the
/// right edge as built by SlabBuilder, a positive value heats the slab (verified
/// against the series below; the opposite of the left-edge Topaz2D_FluxBC golden,
/// where a negative value heats). Analytical solution:
/// Carslaw & Jaeger p. 112. Expected values are the exact series at the report's
/// checkpoints (hygrothermfem_python, analytic.slab_constant_flux); the tolerance
/// is the measured backward-Euler discretization band at dt = 0.01 on ten
/// elements (max 6.8e-3).
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_FluxBC, TransientSeries)
{
    SCOPED_TRACE("Begin Test: unit slab, constant surface flux.");

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

    // On this segment a positive flux heats the slab -- the report's q = 1.
    constexpr auto surfaceFlux = 1.0;
    multiDomain.thermal().createBC_FixedFlux(21, 22, surfaceFlux);

    constexpr auto dTime = 0.01;
    constexpr auto nSteps = 200;

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

    TestHelper::CsvDump dump("topaz_transient_flux.csv", 11);
    for(std::size_t step = 0; step < solution.size(); ++step)
    {
        dump.addRow(step + 1, TestHelper::bottomRow(solution[step], 11, 2));
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.007885293, 0.059310894, 0.356826246}},
      {50, {0.334790713, 0.458333333, 0.831875953}},
      {100, {0.833343815, 0.958333333, 1.333322852}},
      {200, {1.833333334, 1.958333333, 2.333333333}},
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
