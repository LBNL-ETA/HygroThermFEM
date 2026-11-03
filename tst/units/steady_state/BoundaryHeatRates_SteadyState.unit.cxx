#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"

// Same two-element model as ConvectionBC_2D_SteadyStateThermalDomain: convection on the two
// vertical faces, natural (adiabatic) horizontal faces. The consistent boundary heat rate of a
// convection element with equal nodal temperatures reduces to h * L * (Tsurface - Tair), heat
// leaving the domain positive, and the two faces must balance exactly at convergence.
TEST(BoundaryHeatRates_SteadyState, ConvectionBothSides)
{
    SCOPED_TRACE("Begin Test: Boundary heat rates from a steady-state solution.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    multiDomain.nodes().createNode({.index = 1, .x = 15, .y = 5});
    multiDomain.nodes().createNode({.index = 2, .x = 15, .y = 0});
    multiDomain.nodes().createNode({.index = 3, .x = 5, .y = 5});
    multiDomain.nodes().createNode({.index = 4, .x = 5, .y = 0});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());

    multiDomain.createElement(
      {.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement(
      {.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    constexpr auto hc1 = 20.0;
    constexpr auto temperatureAir1 = -18.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1};

    constexpr auto hc2 = 2.4;
    constexpr auto temperatureAir2 = 21.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2};

    multiDomain.thermal().createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.thermal().createBC_FixedHc(6, 5, bcCoeff2);

    // MultiDomain::steadyState (THERM's call path) pushes the solved temperatures into the node
    // pool; boundaryHeatRates reads the node properties, so the bare IDomain::steadyState (which
    // only returns the vector) is not enough here.
    const auto solution = multiDomain.steadyState().temperature;

    const auto heatRates =
      multiDomain.thermal().boundaryHeatRates(HygroThermFEM::Variable::temperature);

    ASSERT_EQ(heatRates.size(), 2u);

    // Right face (nodes 1-2): surface warmer than the -18 C air, heat leaves the domain.
    // Both nodes solve to the same temperature, so Q = h * L * (Tsurface - Tair).
    constexpr auto edgeLength = 5.0;
    const auto expectedRight = hc1 * edgeLength * (solution[0] - temperatureAir1);
    EXPECT_EQ(heatRates[0].node1, 1u);
    EXPECT_EQ(heatRates[0].node2, 2u);
    EXPECT_NEAR(expectedRight, heatRates[0].heatRate, 1e-9);
    EXPECT_NEAR(12.607759, heatRates[0].heatRate, 1e-4);

    // Left face (nodes 6-5): surface cooler than the 21 C air, heat enters (negative rate).
    const auto expectedLeft = hc2 * edgeLength * (solution[4] - temperatureAir2);
    EXPECT_EQ(heatRates[1].node1, 6u);
    EXPECT_EQ(heatRates[1].node2, 5u);
    EXPECT_NEAR(expectedLeft, heatRates[1].heatRate, 1e-9);
    EXPECT_NEAR(-12.607759, heatRates[1].heatRate, 1e-4);

    // Discrete energy balance: everything that enters leaves.
    EXPECT_NEAR(0.0, heatRates[0].heatRate + heatRates[1].heatRate, 1e-9);
}
