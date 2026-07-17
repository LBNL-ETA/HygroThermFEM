#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::State;

// currentStateSolution() is the "timestep 0" frame of a transient series: a zero-duration
// snapshot of the node state taken before any step, with identically zero fluxes.
TEST(CurrentStateSolution, CapturesInitialConditions)
{
    constexpr double initialTemperature = 21.0;
    constexpr double initialHumidity = 0.5;

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

    const State initialState{.temperature = initialTemperature,
                             .humidity = initialHumidity,
                             .pressure = 101325.0,
                             .liquidPercent = 1.0};

    multiDomain.nodes().createNode({.index = 1, .x = 1, .y = 5, .state = initialState});
    multiDomain.nodes().createNode({.index = 2, .x = 1, .y = 0, .state = initialState});
    multiDomain.nodes().createNode({.index = 3, .x = 0.5, .y = 5, .state = initialState});
    multiDomain.nodes().createNode({.index = 4, .x = 0.5, .y = 0, .state = initialState});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5, .state = initialState});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = initialState});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::TestMaterial());
    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    const auto snapshot = multiDomain.currentStateSolution();

    EXPECT_EQ(snapshot.dTime, 0.0);
    ASSERT_EQ(snapshot.temperature.size(), 6u);
    ASSERT_EQ(snapshot.humidity.size(), 6u);
    ASSERT_EQ(snapshot.waterContent.size(), 6u);
    ASSERT_EQ(snapshot.heatFlux.size(), 6u);
    ASSERT_EQ(snapshot.waterFlux.size(), 6u);

    // The water content of the initial state, evaluated by the engine's own material model:
    // the snapshot must agree with the node-state view exactly.
    const auto expectedWater = multiDomain.nodes().properties(HygroThermFEM::Variable::water);

    for(std::size_t index = 0u; index < 6u; ++index)
    {
        EXPECT_EQ(snapshot.temperature[index], initialTemperature) << "node " << index;
        EXPECT_EQ(snapshot.humidity[index], initialHumidity) << "node " << index;
        EXPECT_EQ(snapshot.waterContent[index], expectedWater[index]) << "node " << index;
        EXPECT_EQ(snapshot.heatFlux[index].x, 0.0);
        EXPECT_EQ(snapshot.heatFlux[index].y, 0.0);
        EXPECT_EQ(snapshot.waterFlux[index].x, 0.0);
        EXPECT_EQ(snapshot.waterFlux[index].y, 0.0);
    }
}
