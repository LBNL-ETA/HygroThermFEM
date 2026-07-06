#include <cmath>
#include <map>
#include <vector>
#include <gtest/gtest.h>

#include "BoundaryCondition2DThermal.hxx"
#include "EnclosureRadiation.hxx"
#include "Nodes.hxx"


using namespace HygroThermFEM;

// Standalone tests for the enclosure radiosity coordinator (B1). It consumes the WCE view-factor
// engine and returns, per segment, the effective radiant temperature of its surroundings.

namespace
{
    // Four corners of the unit square (indices 1..4), winding (0,0)->(0,1)->(1,1)->(1,0) so the
    // segment normals point inward, matching the known unit-square view factors.
    Nodes unitSquareNodes(double temperature1, double temperature2, double temperature3,
                          double temperature4)
    {
        Nodes nodes;
        nodes.createNode({.index = 1, .x = 0.0, .y = 0.0, .state = State{.temperature = temperature1}});
        nodes.createNode({.index = 2, .x = 0.0, .y = 1.0, .state = State{.temperature = temperature2}});
        nodes.createNode({.index = 3, .x = 1.0, .y = 1.0, .state = State{.temperature = temperature3}});
        nodes.createNode({.index = 4, .x = 1.0, .y = 0.0, .state = State{.temperature = temperature4}});
        return nodes;
    }

    std::vector<EnclosureRadiationSegment> closedSquareSegments(double emissivity)
    {
        return {{.node1 = 1, .node2 = 2, .emissivity = emissivity, .enclosureId = 0},
                {.node1 = 2, .node2 = 3, .emissivity = emissivity, .enclosureId = 0},
                {.node1 = 3, .node2 = 4, .emissivity = emissivity, .enclosureId = 0},
                {.node1 = 4, .node2 = 1, .emissivity = emissivity, .enclosureId = 0}};
    }
}   // namespace

// A closed isothermal enclosure: every segment sees its own temperature, so the effective radiant
// temperature equals the surface temperature (and the net radiative flux would be zero).
TEST(EnclosureRadiation, IsothermalClosedSquareSeesSurfaceTemperature)
{
    constexpr double temperature = 30.0;
    auto nodes = unitSquareNodes(temperature, temperature, temperature, temperature);

    EnclosureRadiation enclosure(nodes, closedSquareSegments(0.9), {});

    ASSERT_EQ(4u, enclosure.numberOfSegments());
    for(std::size_t idx = 0; idx < enclosure.numberOfSegments(); ++idx)
    {
        EXPECT_NEAR(temperature, enclosure.effectiveRadiantTemperature(idx), 1e-6);
    }
}

// The two surface-temperature conventions: with a temperature gradient along a segment, the
// Conrad-compatible SegmentIsothermal model radiates at the fourth-power mean of the node
// temperatures (evaluated in Kelvin), while LocalTemperature uses the linear mean. Equal node
// temperatures make the conventions coincide.
TEST(EnclosureRadiation, SurfaceTemperatureConventions)
{
    constexpr double coldEnd = 0.0;
    constexpr double hotEnd = 20.0;
    auto nodes = unitSquareNodes(coldEnd, hotEnd, hotEnd, coldEnd);
    const auto segments = closedSquareSegments(0.9);

    EnclosureRadiation isothermal(
      nodes, segments, {}, true, EnclosureSurfaceTemperature::SegmentIsothermal);
    EnclosureRadiation local(
      nodes, segments, {}, true, EnclosureSurfaceTemperature::LocalTemperature);

    // Left wall (nodes 1-2) spans 0 C to 20 C.
    constexpr double kelvinOffset = 273.15;
    const double fourthPowerMean =
      std::pow(0.5 * (std::pow(coldEnd + kelvinOffset, 4) + std::pow(hotEnd + kelvinOffset, 4)),
               0.25)
      - kelvinOffset;

    EXPECT_NEAR(fourthPowerMean, isothermal.segmentSurfaceTemperature(0), 1e-9);
    EXPECT_NEAR(0.5 * (coldEnd + hotEnd), local.segmentSurfaceTemperature(0), 1e-9);
    EXPECT_GT(isothermal.segmentSurfaceTemperature(0), local.segmentSurfaceTemperature(0));

    // Top wall (nodes 2-3) is uniform at 20 C: the conventions coincide.
    EXPECT_NEAR(hotEnd, isothermal.segmentSurfaceTemperature(1), 1e-9);
    EXPECT_NEAR(hotEnd, local.segmentSurfaceTemperature(1), 1e-9);
}

// A closed enclosure with hot (left) and cold (right) walls: each segment's radiant temperature
// lies between the extremes, and a wall sees surroundings colder than itself if it is hot and
// warmer if it is cold (so heat flows hot -> cold).
TEST(EnclosureRadiation, TwoTemperatureClosedSquareFlowsHotToCold)
{
    constexpr double hot = 40.0;
    constexpr double cold = 20.0;
    // Left wall (nodes 1,2) hot; right wall (nodes 3,4) cold; top/bottom averages 30.
    auto nodes = unitSquareNodes(hot, hot, cold, cold);

    EnclosureRadiation enclosure(nodes, closedSquareSegments(0.9), {});

    const auto leftWall = enclosure.effectiveRadiantTemperature(0);    // nodes 1-2, hot
    const auto rightWall = enclosure.effectiveRadiantTemperature(2);   // nodes 3-4, cold

    EXPECT_GT(leftWall, cold);
    EXPECT_LT(leftWall, hot);     // hot wall sees colder surroundings -> loses heat
    EXPECT_LT(rightWall, hot);
    EXPECT_GT(rightWall, cold);   // cold wall sees warmer surroundings -> gains heat
}

// An open (auto) enclosure: three walls open to an environment. When the surfaces and the
// environment are all at the same temperature, every segment sees that temperature.
TEST(EnclosureRadiation, IsothermalOpenEnclosureSeesEnvironmentTemperature)
{
    constexpr double temperature = 25.0;
    auto nodes = unitSquareNodes(temperature, temperature, temperature, temperature);

    // Three walls of the square (left, right, bottom), open at the top; enclosure id 1.
    const std::vector<EnclosureRadiationSegment> walls{
      {.node1 = 1, .node2 = 2, .emissivity = 0.84, .enclosureId = 1},   // left
      {.node1 = 3, .node2 = 4, .emissivity = 0.84, .enclosureId = 1},   // right
      {.node1 = 4, .node2 = 1, .emissivity = 0.84, .enclosureId = 1}};  // bottom

    EnclosureRadiation enclosure(nodes, walls, {{1u, temperature}});

    for(std::size_t idx = 0; idx < enclosure.numberOfSegments(); ++idx)
    {
        EXPECT_NEAR(temperature, enclosure.effectiveRadiantTemperature(idx), 1e-6);
    }
}

// EnclosureRadiationBC: at an isothermal enclosure each segment's radiant temperature equals its
// surface temperature, so the linearized contribution q = h (T_rad - T_s) is zero. In assembly
// terms the right-hand side equals the H-matrix times the surface temperatures.
TEST(EnclosureRadiationBC, IsothermalEnclosureGivesZeroNetFlux)
{
    constexpr double temperature = 30.0;
    constexpr double emissivity = 0.9;
    auto nodes = unitSquareNodes(temperature, temperature, temperature, temperature);

    EnclosureRadiation coordinator(nodes, closedSquareSegments(emissivity), {});
    EnclosureRadiationBC boundary(nodes, 1, 2, emissivity, coordinator, 0);

    const auto rightHandSide = boundary.R_Vector();
    const auto hMatrix = boundary.H_Matrix();

    const std::vector<double> surfaceTemperatures{temperature, temperature};
    for(std::size_t row = 0; row < 2; ++row)
    {
        double hTimesTemperature = 0.0;
        for(std::size_t col = 0; col < 2; ++col)
        {
            hTimesTemperature += hMatrix(row, col) * surfaceTemperatures[col];
        }
        EXPECT_NEAR(rightHandSide[row], hTimesTemperature, 1e-3);
    }
}
