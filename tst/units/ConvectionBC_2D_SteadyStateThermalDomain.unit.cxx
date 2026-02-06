#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

TEST(ConvectionBC_2D_SteadyStateThermalDomain, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with simple conduction.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    multiDomain.nodes().createNode({.index = 1, .x = 15, .y = 5});
    multiDomain.nodes().createNode({.index = 2, .x = 15, .y = 0});
    multiDomain.nodes().createNode({.index = 3, .x = 5, .y = 5});
    multiDomain.nodes().createNode({.index = 4, .x = 5, .y = 0});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto hc1 = 20.0;
    constexpr auto temperatureAir1 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1};

    constexpr auto hc2 = 2.4;
    constexpr auto temperatureAir2 = 21.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2};

    multiDomain.thermal().createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.thermal().createBC_FixedHc(6, 5, bcCoeff2);

    auto solution = multiDomain.thermal().steadyState();

    std::vector<double> correctSolution{
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        EXPECT_NEAR(correctSolution[i], solution[i], 1e-6);
    }
}
