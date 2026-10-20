#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "EnclosureRadiation.hxx"
#include "TestMaterials.hxx"

using namespace HygroThermFEM;

// The steady pseudo-time march must stop only when the field satisfies the nonlinear steady
// equations. The model is the two-bar radiating gap of EnclosureRadiation_2D_SteadyState: the
// bars are coupled only by enclosure radiation, so every pass reassembles a temperature-dependent
// system.
namespace
{
    constexpr double tolerance = 1e-8;

    // MultiDomain's domains reference its own node and material pools, so it is built in place
    // rather than returned by value.
    void buildTwoBarsRadiatingGap(MultiDomain & multiDomain, const size_t maxIterations)
    {
        multiDomain.setSolverSettings({.relaxationParameter = 1.0,
                                       .errorTolerance = tolerance,
                                       .maxNumberOfIterations = maxIterations,
                                       .maxDivisions = 3,
                                       .numberOfSubtimesteps = 10});

        auto & nodes = multiDomain.nodes();
        nodes.createNode({.x = 0.0, .y = 0.0});   // 1
        nodes.createNode({.x = 1.0, .y = 0.0});   // 2  bar 1 right face (with 3)
        nodes.createNode({.x = 1.0, .y = 1.0});   // 3
        nodes.createNode({.x = 0.0, .y = 1.0});   // 4
        nodes.createNode({.x = 2.0, .y = 0.0});   // 5  bar 2 left face (with 8)
        nodes.createNode({.x = 3.0, .y = 0.0});   // 6
        nodes.createNode({.x = 3.0, .y = 1.0});   // 7
        nodes.createNode({.x = 2.0, .y = 1.0});   // 8

        const auto & material =
          multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());
        multiDomain.createElement(
          {.node1 = 1, .node2 = 2, .node3 = 3, .node4 = 4, .material = material.name()});
        multiDomain.createElement(
          {.node1 = 5, .node2 = 6, .node3 = 7, .node4 = 8, .material = material.name()});

        multiDomain.createBC_FixedTemperature(1, 4, 40.0);
        multiDomain.createBC_FixedTemperature(6, 7, 20.0);

        multiDomain.createEnclosureRadiation(
          {{.node1 = 2, .node2 = 3, .emissivity = 0.9, .enclosureId = 0},
           {.node1 = 8, .node2 = 5, .emissivity = 0.9, .enclosureId = 0}},
          {{0u, 30.0}});
    }
}   // namespace

TEST(SteadyStateConvergenceCheck, ConvergedFieldSatisfiesSteadyEquations)
{
    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    buildTwoBarsRadiatingGap(multiDomain, 200);

    const auto solution = multiDomain.steadyState();

    // The pseudo-time step of the march's steady regime.
    constexpr double steadyPseudoStep = 1.0e12;
    EXPECT_TRUE(solution.converged);
    EXPECT_LE(solution.temperatureError, tolerance);
    EXPECT_LE(multiDomain.thermal().steadyStateResidual(steadyPseudoStep), tolerance);
}
