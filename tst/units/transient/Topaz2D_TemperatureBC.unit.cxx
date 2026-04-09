#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 100 degrees
///   Initial temperature boundary conditions at nodes 5 and 6 are 12 degrees
///   Solution achieved with linear solver (no iterations required in this case
/////////////////////////////////////////////////////////////////////////////////////

TEST(Topaz2D_TemperatureBC, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with transient.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate, initial temperature

    // same temperature in every node (humidity and pressure irrelevant for this example)
    const auto state = HygroThermFEM::State{
        .temperature = 100,
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
    constexpr auto tSurface = 12.0;

    multiDomain.thermal().createBC_FixedTemperature(5, 6, tSurface);

    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, 3600, 4);

    std::vector<std::vector<double>> correctSolution{
      {83.64609365, 83.64609365, 61.65791323, 61.65791323, 12, 12},
      {66.21082587, 66.21082587, 42.76873166, 42.76873166, 12, 12},
      {51.74326318, 51.74326318, 32.29131256, 32.29131256, 12, 12},
      {40.71210006, 40.71210006, 25.88046294, 25.88046294, 12, 12}};

    TestHelper::expectNear(correctSolution, solution, 1e-6);
}
