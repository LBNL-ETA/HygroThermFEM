#include <array>
#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.4a -- Plane Transient Linear, Constant Element Generation.
///
/// Unit slab (alpha = 1, k = 1) at zero initial temperature, adiabatic at x = 0,
/// surface at x = 1 held at 0, uniform generation q = 1 switched on at t = 0.
/// Analytical solution: the steady parabola minus its decaying odd-cosine series
/// (report section 3.4). Expected values are the exact series at the report's
/// checkpoints (hygrothermfem_python, analytic.slab_generation_transient); the
/// tolerance is the measured backward-Euler discretization band at dt = 0.01 on
/// ten elements (max 2.7e-3).
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_GenerationTransient, UniformGeneration)
{
    SCOPED_TRACE("Begin Test: unit slab, uniform generation switched on at t = 0.");

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

    constexpr auto tSurface = 0.0;
    multiDomain.thermal().createBC_FixedTemperature(21, 22, tSurface);

    constexpr auto generation = 1.0;
    multiDomain.thermal().setVolumetricSource(generation);

    constexpr auto dTime = 0.01;
    constexpr auto nSteps = 200;

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

    TestHelper::CsvDump dump("topaz_generation_transient.csv", 11);
    for(std::size_t step = 0; step < solution.size(); ++step)
    {
        dump.addRow(step + 1, TestHelper::bottomRow(solution[step], 11, 2));
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.098873183, 0.088439135, 0.0}},
      {50, {0.349727265, 0.268740723, 0.0}},
      {100, {0.456238552, 0.344055983, 0.0}},
      {200, {0.496288812, 0.372375794, 0.0}},
    };

    for(const auto & [step, expected] : checkpoints)
    {
        for(std::size_t pos = 0; pos < expected.size(); ++pos)
        {
            EXPECT_NEAR(expected[pos], solution[step - 1][pos * 10], 0.005)
              << "step " << step << ", position " << pos;
        }
    }
}
