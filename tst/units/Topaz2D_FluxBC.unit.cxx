#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

TEST(Topaz2D_FluxBC, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    const auto state = HygroThermFEM::State{
        .temperature = 0,
        .humidity = 0,
        .pressure = 101325,
        .liquidPercent = 0
    };

    multiDomain.nodes().createNode({.index = 1, .x = 0.15, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 0.15, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 0.05, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 0.05, .y = 0, .state = state});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 0.05, .state = state});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = state});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::CottaerSandstoneNonPorous());

    multiDomain.createElement({.node1 = 1, .node2 = 2, .node3 = 4, .node4 = 3, .material = material.name()});
    multiDomain.createElement({.node1 = 5, .node2 = 3, .node3 = 4, .node4 = 6, .material = material.name()});

    // Create Boundary Conditions
    // Positive flux means outside flow.
    constexpr auto surfaceFlux = -12.0;

    multiDomain.thermal().createBC_FixedFlux(5, 6, surfaceFlux);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;

    auto temperatures = multiDomain.nodes().properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = multiDomain.thermal().transient(temperatures, dTime).solution;
        solution.push_back(temperatures);
    }

    std::vector<std::vector<double>> correctSolution{
      {0.068797095, 0.068797095, 0.161296275, 0.161296275, 0.370195609, 0.370195609},
      {0.184225668, 0.184225668, 0.339421878, 0.339421878, 0.596640275, 0.596640275},
      {0.324790684, 0.324790684, 0.513783385, 0.513783385, 0.784104345, 0.784104345},
      {0.478007410, 0.478007410, 0.684010609, 0.684010609, 0.958667844, 0.958667844}};

    TestHelper::expectNear(correctSolution, solution, 1e-6);
}
