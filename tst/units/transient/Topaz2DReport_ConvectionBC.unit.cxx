#include <array>
#include <gtest/gtest.h>

#include "DumpCsv.hxx"
#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// TOPAZ2D Validation Status Report (Drayton, UCRL-ID-106550, December 1990),
/// section 3.3 -- Plane Transient Linear, Constant Convection Boundary Condition.
///
/// Unit slab (alpha = 1, k = 1) at uniform initial temperature 1.0, adiabatic at
/// x = 0, convection h (T - 0) with h = 1 at x = 1. Analytical solution: Biot
/// series with the roots of beta tan(beta) = h L / k, Carslaw & Jaeger p. 122.
/// Expected values are the exact series at the report's checkpoints
/// (hygrothermfem_python, analytic.slab_convection); the tolerance is the
/// measured backward-Euler discretization band at dt = 0.01 on ten elements
/// (max 4.7e-3).
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2DReport_ConvectionBC, TransientSeries)
{
    SCOPED_TRACE("Begin Test: unit slab, constant convection surface.");

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

    constexpr auto tAir = 0.0;
    constexpr auto hc = 1.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tAir, hc};
    multiDomain.thermal().createBC_FixedHc(21, 22, bcCoeff);

    constexpr auto dTime = 0.01;
    constexpr auto nSteps = 200;

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

    TestHelper::CsvDump dump("topaz_transient_convection.csv", 11);
    for(std::size_t step = 0; step < solution.size(); ++step)
    {
        dump.addRow(step + 1, TestHelper::bottomRow(solution[step], 11, 2));
    }

    // Exact series at x = 0, 0.5, 1.0 for steps 10, 50, 100, 200 (t = step * dt).
    const std::vector<std::pair<std::size_t, std::array<double, 3>>> checkpoints{
      {10, {0.993108255, 0.950508452, 0.723577239}},
      {50, {0.772526383, 0.702597259, 0.504521928}},
      {100, {0.533859401, 0.485224060, 0.348176852}},
      {200, {0.254668042, 0.231466817, 0.166090581}},
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
