#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "EnclosureRadiation.hxx"
#include "TestMaterials.hxx"

using namespace HygroThermFEM;

namespace
{
    //! Two conducting bars separated by a gap, as in EnclosureRadiation_2D_SteadyState: the
    //! far faces are held at 40 C and 20 C, and the facing surfaces form an open enclosure
    //! that radiates to each other and to a 30 C environment. Radiation is the only
    //! coupling between the bars.
    void buildTwoBars(MultiDomain & multiDomain)
    {
        multiDomain.nodes().createNode({.x = 0.0, .y = 0.0});   // 1
        multiDomain.nodes().createNode({.x = 1.0, .y = 0.0});   // 2  bar 1 right face (with 3)
        multiDomain.nodes().createNode({.x = 1.0, .y = 1.0});   // 3
        multiDomain.nodes().createNode({.x = 0.0, .y = 1.0});   // 4
        multiDomain.nodes().createNode({.x = 2.0, .y = 0.0});   // 5  bar 2 left face (with 8)
        multiDomain.nodes().createNode({.x = 3.0, .y = 0.0});   // 6
        multiDomain.nodes().createNode({.x = 3.0, .y = 1.0});   // 7
        multiDomain.nodes().createNode({.x = 2.0, .y = 1.0});   // 8

        const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());
        multiDomain.createElement({.node1 = 1, .node2 = 2, .node3 = 3, .node4 = 4, .material = material.name()});
        multiDomain.createElement({.node1 = 5, .node2 = 6, .node3 = 7, .node4 = 8, .material = material.name()});

        multiDomain.createBC_FixedTemperature(1, 4, 40.0);
        multiDomain.createBC_FixedTemperature(6, 7, 20.0);

        multiDomain.createEnclosureRadiation({{.node1 = 2, .node2 = 3, .emissivity = 0.9, .enclosureId = 0},
                                              {.node1 = 8, .node2 = 5, .emissivity = 0.9, .enclosureId = 0}},
                                             {{0u, 30.0}});
    }
}   // namespace

// An enclosure radiates in a transient run as it does in a steady one: marched to rest
// under the same constant conditions, the transient reaches the steady solution.
TEST(EnclosureRadiation_2D_Transient, TwoBarsReachTheSteadySolution)
{
    MultiDomain steadyDomain({.performThermal = true, .performMoisture = false});
    buildTwoBars(steadyDomain);
    const auto steady = steadyDomain.steadyState();
    ASSERT_TRUE(steady.converged);

    MultiDomain transientDomain({.performThermal = true, .performMoisture = false});
    buildTwoBars(transientDomain);
    auto temperatures = transientDomain.nodes().properties(Variable::temperature);
    auto humidities = transientDomain.nodes().properties(Variable::humidity);
    constexpr double stepSeconds{1.0e6};
    constexpr std::size_t steps{400u};
    for(std::size_t step = 0u; step < steps; ++step)
    {
        const auto solution = transientDomain.transient(temperatures, humidities, stepSeconds, step);
        temperatures = solution.temperature;
        humidities = solution.humidity;
    }

    ASSERT_EQ(steady.temperature.size(), temperatures.size());
    for(std::size_t node = 0u; node < temperatures.size(); ++node)
    {
        EXPECT_NEAR(temperatures[node], steady.temperature[node], 1e-3) << "node " << node + 1u;
    }

    // The radiating faces did move: without the enclosure each bar would sit at its own
    // fixed temperature throughout.
    EXPECT_LT(0.5 * (temperatures[1] + temperatures[2]), 39.9);
    EXPECT_GT(0.5 * (temperatures[4] + temperatures[7]), 20.1);
}
