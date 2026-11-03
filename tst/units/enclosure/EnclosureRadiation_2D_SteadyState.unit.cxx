#include <cmath>
#include <map>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "EnclosureRadiation.hxx"
#include "TestMaterials.hxx"


using namespace HygroThermFEM;

// B3: end-to-end steady-state with enclosure radiation. Two conducting bars separated by a gap;
// their facing surfaces form an open (auto) enclosure that radiates to each other and to a 30 C
// environment. Bar 1's far side is held at 40 C, bar 2's at 20 C. The bars are thermally
// disconnected by conduction (the gap carries no elements), so the only coupling is the nonlinear
// radiation. The steady solver drives the (non-linear) enclosure BCs to convergence.
TEST(EnclosureRadiation_2D_SteadyState, TwoBarsRadiatingGap)
{
    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Bar 1: [0,1] x [0,1] (nodes 1..4); Bar 2: [2,3] x [0,1] (nodes 5..8).
    multiDomain.nodes().createNode({.x = 0.0, .y = 0.0});   // 1
    multiDomain.nodes().createNode({.x = 1.0, .y = 0.0});   // 2  bar 1 right face (with 3)
    multiDomain.nodes().createNode({.x = 1.0, .y = 1.0});   // 3
    multiDomain.nodes().createNode({.x = 0.0, .y = 1.0});   // 4
    multiDomain.nodes().createNode({.x = 2.0, .y = 0.0});   // 5  bar 2 left face (with 8)
    multiDomain.nodes().createNode({.x = 3.0, .y = 0.0});   // 6
    multiDomain.nodes().createNode({.x = 3.0, .y = 1.0});   // 7
    multiDomain.nodes().createNode({.x = 2.0, .y = 1.0});   // 8

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());
    multiDomain.createElement(
      {.node1 = 1, .node2 = 2, .node3 = 3, .node4 = 4, .material = material.name()});
    multiDomain.createElement(
      {.node1 = 5, .node2 = 6, .node3 = 7, .node4 = 8, .material = material.name()});

    // Fixed temperatures on the outer faces.
    multiDomain.createBC_FixedTemperature(1, 4, 40.0);   // bar 1 left
    multiDomain.createBC_FixedTemperature(6, 7, 20.0);   // bar 2 right

    // Open enclosure (id 0) between the facing surfaces, environment at 30 C.
    // Segment A = bar 1 right face (2->3, normal +x); segment B = bar 2 left face (8->5, normal -x).
    multiDomain.createEnclosureRadiation(
      {{.node1 = 2, .node2 = 3, .emissivity = 0.9, .enclosureId = 0},
       {.node1 = 8, .node2 = 5, .emissivity = 0.9, .enclosureId = 0}},
      {{0u, 30.0}});

    const auto solution = multiDomain.steadyState();
    const auto & temperature = solution.temperature;

    ASSERT_EQ(8u, temperature.size());
    for(const auto value : temperature)
    {
        EXPECT_TRUE(std::isfinite(value));   // converged, did not diverge
    }

    // Fixed faces hold their setpoints.
    EXPECT_NEAR(40.0, temperature[0], 0.1);   // node 1
    EXPECT_NEAR(40.0, temperature[3], 0.1);   // node 4
    EXPECT_NEAR(20.0, temperature[5], 0.1);   // node 6
    EXPECT_NEAR(20.0, temperature[6], 0.1);   // node 7

    // The facing faces are coupled only by radiation: both lie strictly between the extremes, and
    // bar 1's face stays hotter than bar 2's (heat flows bar 1 -> bar 2).
    const double bar1Face = 0.5 * (temperature[1] + temperature[2]);   // nodes 2, 3
    const double bar2Face = 0.5 * (temperature[4] + temperature[7]);   // nodes 5, 8
    EXPECT_GT(bar1Face, 20.0);
    EXPECT_LT(bar1Face, 40.0);
    EXPECT_GT(bar2Face, 20.0);
    EXPECT_LT(bar2Face, 40.0);
    EXPECT_GT(bar1Face, bar2Face);
}

// Closed frame cavity: four conducting bars surround a central [0,1]x[0,1] cavity, and their inner
// faces form a CLOSED radiation enclosure (the unit-square view factors). Each bar's outer face is
// held at a different temperature (bottom 40, top 20, sides 30). The bars are conduction-isolated
// from each other; the closed-enclosure radiosity redistributes heat across the cavity walls. The
// steady solver converges and the cavity walls settle between the extremes, hottest at the bottom.
TEST(EnclosureRadiation_2D_SteadyState, ClosedFrameCavity)
{
    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    auto & nodes = multiDomain.nodes();
    // Bottom bar (nodes 1..4): x[0,1] y[-0.5,0], inner (top) face = cavity bottom wall.
    nodes.createNode({.x = 0.0, .y = -0.5});   // 1
    nodes.createNode({.x = 1.0, .y = -0.5});   // 2
    nodes.createNode({.x = 1.0, .y = 0.0});    // 3
    nodes.createNode({.x = 0.0, .y = 0.0});    // 4
    // Top bar (5..8): x[0,1] y[1,1.5], inner (bottom) face = cavity top wall.
    nodes.createNode({.x = 0.0, .y = 1.0});    // 5
    nodes.createNode({.x = 1.0, .y = 1.0});    // 6
    nodes.createNode({.x = 1.0, .y = 1.5});    // 7
    nodes.createNode({.x = 0.0, .y = 1.5});    // 8
    // Left bar (9..12): x[-0.5,0] y[0,1], inner (right) face = cavity left wall.
    nodes.createNode({.x = -0.5, .y = 0.0});   // 9
    nodes.createNode({.x = 0.0, .y = 0.0});    // 10
    nodes.createNode({.x = 0.0, .y = 1.0});    // 11
    nodes.createNode({.x = -0.5, .y = 1.0});   // 12
    // Right bar (13..16): x[1,1.5] y[0,1], inner (left) face = cavity right wall.
    nodes.createNode({.x = 1.0, .y = 0.0});    // 13
    nodes.createNode({.x = 1.5, .y = 0.0});    // 14
    nodes.createNode({.x = 1.5, .y = 1.0});    // 15
    nodes.createNode({.x = 1.0, .y = 1.0});    // 16

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());
    multiDomain.createElement({.node1 = 1, .node2 = 2, .node3 = 3, .node4 = 4, .material = material.name()});
    multiDomain.createElement({.node1 = 5, .node2 = 6, .node3 = 7, .node4 = 8, .material = material.name()});
    multiDomain.createElement({.node1 = 9, .node2 = 10, .node3 = 11, .node4 = 12, .material = material.name()});
    multiDomain.createElement({.node1 = 13, .node2 = 14, .node3 = 15, .node4 = 16, .material = material.name()});

    // Outer faces held at fixed temperatures.
    multiDomain.createBC_FixedTemperature(1, 2, 40.0);     // bottom bar outer
    multiDomain.createBC_FixedTemperature(7, 8, 20.0);     // top bar outer
    multiDomain.createBC_FixedTemperature(9, 12, 30.0);    // left bar outer
    multiDomain.createBC_FixedTemperature(14, 15, 30.0);   // right bar outer

    // Closed enclosure (id 0): the four cavity walls, each with an inward normal.
    multiDomain.createEnclosureRadiation(
      {{.node1 = 3, .node2 = 4, .emissivity = 0.9, .enclosureId = 0},     // bottom wall
       {.node1 = 10, .node2 = 11, .emissivity = 0.9, .enclosureId = 0},   // left wall
       {.node1 = 5, .node2 = 6, .emissivity = 0.9, .enclosureId = 0},     // top wall
       {.node1 = 16, .node2 = 13, .emissivity = 0.9, .enclosureId = 0}}); // right wall

    const auto solution = multiDomain.steadyState();
    const auto & temperature = solution.temperature;

    ASSERT_EQ(16u, temperature.size());
    for(const auto value : temperature)
    {
        EXPECT_TRUE(std::isfinite(value));
    }

    const double bottomWall = 0.5 * (temperature[2] + temperature[3]);    // nodes 3,4
    const double topWall = 0.5 * (temperature[4] + temperature[5]);       // nodes 5,6
    const double leftWall = 0.5 * (temperature[9] + temperature[10]);     // nodes 10,11
    const double rightWall = 0.5 * (temperature[15] + temperature[12]);   // nodes 16,13

    for(const double wall : {bottomWall, topWall, leftWall, rightWall})
    {
        EXPECT_GT(wall, 20.0);
        EXPECT_LT(wall, 40.0);
    }
    EXPECT_GT(bottomWall, topWall);          // bottom (40 C side) hotter than top (20 C side)
    EXPECT_LT(bottomWall, 40.0);             // radiation cools the hot wall
    EXPECT_GT(topWall, 20.0);                // radiation warms the cold wall
}
