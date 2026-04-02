#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(ConvectionBC_2D_TransientNoChanges, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    constexpr HygroThermFEM::State state({
        .temperature = 20.0,
        .humidity = 0.0,
        .pressure = 101325.0,
        .liquidPercent = 0
    });

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstone());

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1, 0.15})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    // Create Boundary Conditions
    constexpr auto tSurface = 20.0;
    constexpr auto hc = 1.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tSurface, hc};

    multiDomain.thermal().createBC_FixedHc(1, 2, bcCoeff);

    constexpr auto dTime = 36000;
    constexpr auto nSteps = 4;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0},
      {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}};

    TestHelper::expectNear(correctSolution, solution, 1e-6);
}
