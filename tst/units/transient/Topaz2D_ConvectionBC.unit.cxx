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

/// Transient convection BC. Example is tested in EXCEL and Topaz2D.

TEST(Topaz2D_ConvectionBC, TestExample_1)
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

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto tSurface = 20.0;
    constexpr auto hc = 20.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{tSurface, hc};

    multiDomain.thermal().createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 4;


    const auto solution = multiDomain.thermal().transientMultiStep(
      HygroThermFEM::Variable::temperature, dTime, nSteps);

    // Bottom nodes in ascending x (nodes 6, 4, 2 at x = 0, 0.05, 0.15); the node
    // numbering above runs in descending x, so the row is picked out explicitly.
    for(std::size_t step = 0; step < solution.size(); ++step)
    {
    }

    std::vector<std::vector<double>> correctSolution{
      {1.4182108, 1.4182108, 3.3250259, 3.3250259, 7.6313602, 7.6313602},
      {3.4666897, 3.4666897, 6.2209137, 6.2209137, 10.518214, 10.518214},
      {5.6122051, 5.6122051, 8.4968969, 8.4968969, 12.234324, 12.234324},
      {7.6189241, 7.6189241, 10.317001, 10.317001, 13.501417, 13.501417}};

    TestHelper::expectNear(correctSolution, solution, 1e-6);
}
