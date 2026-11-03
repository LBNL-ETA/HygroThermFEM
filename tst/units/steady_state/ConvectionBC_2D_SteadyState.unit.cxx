#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

TEST(ConvectionBC_2D_SteadyState, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with simple conduction.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    multiDomain.nodes().createNode({.x = 15, .y = 5});
    multiDomain.nodes().createNode({.x = 15, .y = 0});
    multiDomain.nodes().createNode({.x = 5, .y = 5});
    multiDomain.nodes().createNode({.x = 5, .y = 0});
    multiDomain.nodes().createNode({.x = 0, .y = 5});
    multiDomain.nodes().createNode({.x = 0, .y = 0});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    constexpr auto hc1 = 20.0;
    constexpr auto humidity1 = 0.0;
    constexpr auto temperatureAir1 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1, humidity1};

    constexpr auto hc2 = 2.4;
    constexpr auto humidity2 = 0.0;
    constexpr auto temperatureAir2 = 21.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2, humidity2};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = multiDomain.steadyState();
    auto temperature = solution.temperature;

    std::vector<double> correctTemperature{
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};

    EXPECT_EQ(temperature.size(), correctTemperature.size());

    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(temperature[i], correctTemperature[i], 1e-6);
    }
}
