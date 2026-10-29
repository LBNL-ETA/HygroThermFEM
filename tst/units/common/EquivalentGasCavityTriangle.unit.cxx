#include <gtest/gtest.h>
#include <vector>

#include "HygroThermFEM2D.hxx"
#include "GasCavities.hxx"

// A triangular frame cavity has three walls, so one of the four sides the engine
// rectangularizes into has no wall at all. Conrad gives such a side the temperature of the
// cavity node farthest in that direction and the mean emissivity of its two neighbouring
// sides; the engine must do the same instead of evaluating the side at absolute zero.
//
// The triangle: a 16.7 mm by 0.6 mm sliver with a warm node at its left tip and the short
// vertical wall on the right. Its walls are Bottom, Right and Top; Left is empty.

namespace
{
    using HygroThermFEM::cavitySegmentKey;

    constexpr double wallEmissivity{0.9};
    constexpr double warmCelsius{10.0};
    constexpr double coldCelsius{0.0};
    constexpr double zeroCelsiusInKelvin{273.15};

    void createTriangleNodes(HygroThermFEM::Nodes & nodes)
    {
        nodes.createNode({.index = 1u, .x = 0.0, .y = 0.0, .state = {.temperature = warmCelsius}});
        nodes.createNode({.index = 2u, .x = 0.0167, .y = 0.0, .state = {.temperature = coldCelsius}});
        nodes.createNode({.index = 3u, .x = 0.0167, .y = 0.0006, .state = {.temperature = coldCelsius}});
    }

    HygroThermFEM::CavitySurfaceEmissivities allWalls()
    {
        return {{cavitySegmentKey(1u, 2u), wallEmissivity},
                {cavitySegmentKey(2u, 3u), wallEmissivity},
                {cavitySegmentKey(3u, 1u), wallEmissivity}};
    }

    HygroThermFEM::GasCavityState evaluate(const std::vector<std::size_t> & loop)
    {
        HygroThermFEM::Nodes nodes;
        createTriangleNodes(nodes);
        HygroThermFEM::Materials materials;
        materials.createGas("cavity", HygroThermFEM::CavityStandard::ISO15099);
        HygroThermFEM::EquivalentGasCavity cavity(nodes, loop, materials.gas("cavity"), allWalls());
        cavity.update();
        return cavity.state();
    }
}   // namespace

TEST(EquivalentGasCavityTriangle, EmptySideTakesConradsFallback)
{
    const auto state{evaluate({1u, 2u, 3u})};

    // Heat flows from the warm left tip to the cold right wall.
    EXPECT_EQ(KeffCavity::ScreenFlow::Right, state.heatFlow);
    // The empty left side takes the temperature of the leftmost node, not absolute zero,
    // and the emissivity of its neighbours.
    EXPECT_NEAR(warmCelsius + zeroCelsiusInKelvin, state.side1.temperature, 1e-9);
    EXPECT_NEAR(wallEmissivity, state.side1.emissivity, 1e-9);
    // The right wall is a real wall.
    EXPECT_NEAR(coldCelsius + zeroCelsiusInKelvin, state.side2.temperature, 1e-9);
    EXPECT_NEAR(wallEmissivity, state.side2.emissivity, 1e-9);
    EXPECT_GT(state.conductivity, 0.0);
}

// The steady-state physics flags turn the temperature dependence off, so the thermal
// element reads the material's dry conductivity. The evaluated keff must land there too, or
// a steady solve would run the cavity at the gas default and the reported state would
// describe a conductivity the field never saw.
TEST(EquivalentGasCavityTriangle, EvaluatedConductivityReachesTheDryConductivity)
{
    HygroThermFEM::Nodes nodes;
    createTriangleNodes(nodes);
    HygroThermFEM::Materials materials;
    materials.createGas("cavity", HygroThermFEM::CavityStandard::ISO15099);
    HygroThermFEM::EquivalentGasCavity cavity(nodes, {1u, 2u, 3u}, materials.gas("cavity"), allWalls());
    cavity.update();

    EXPECT_NEAR(cavity.state().conductivity, materials.gas("cavity").thermalConductivityDry(), 1e-12);
}

TEST(EquivalentGasCavityTriangle, ClockwiseLoopEvaluatesLikeCounterClockwise)
{
    const auto counterClockwise{evaluate({1u, 2u, 3u})};
    const auto clockwise{evaluate({3u, 2u, 1u})};

    EXPECT_EQ(counterClockwise.heatFlow, clockwise.heatFlow);
    EXPECT_NEAR(counterClockwise.side1.temperature, clockwise.side1.temperature, 1e-9);
    EXPECT_NEAR(counterClockwise.side2.temperature, clockwise.side2.temperature, 1e-9);
    EXPECT_NEAR(counterClockwise.conductivity, clockwise.conductivity, 1e-12);
}
