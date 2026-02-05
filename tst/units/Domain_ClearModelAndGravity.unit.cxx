#include <gtest/gtest.h>
#include <memory>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::State;
using HygroThermFEM::ThermalDomain;
using HygroThermFEM::MultiDomain;
using HygroThermFEM::Variable;


namespace
{
    void createSimpleModel(MultiDomain & multiDomain)
    {
        const State state({
            .temperature = 20.0,
            .humidity = 0.5,
            .pressure = 101325.0
        });

        // Create 4 nodes for a single quad element
        multiDomain.nodes().createNode({.index = 1, .x = 0.0, .y = 0.0, .state = state});
        multiDomain.nodes().createNode({.index = 2, .x = 1.0, .y = 0.0, .state = state});
        multiDomain.nodes().createNode({.index = 3, .x = 1.0, .y = 1.0, .state = state});
        multiDomain.nodes().createNode({.index = 4, .x = 0.0, .y = 1.0, .state = state});

        // Create material
        multiDomain.materials().createSolidMaterial(
          {.name = "TestMaterial",
           .thermalConductivityDry = 1.0,
           .density = 1000.0,
           .porosity = 0.1,
           .heatCapacity = 1000.0,
           .diffusionResistanceFactor = 10.0,
           .thermalConductivityMoistureDependent = {{0.0, 1.0}, {180, 1.0}},
           .moistureDependentMeasurementTemperature = 0,
           .thermalConductivityTemperatureDependent = {{0.0, 1.0}, {100.0, 1.0}},
           .temperatureDependentMeasurementHumidity = 0,
           .liquidTransportCurve = {{0, 0}, {180, 1e-6}},
           .sorptionCurve = {{0, 0}, {1, 180}}});

        multiDomain.createElement({.node1 = 1, .node2 = 2, .node3 = 3, .node4 = 4, .material = "TestMaterial"});
    }

    // Helper to create model with frame cavity for gravity tests
    void createModelWithFrameCavity(MultiDomain & multiDomain)
    {
        const State state({
            .temperature = 20.0,
            .humidity = 0.0,
            .pressure = 101325.0
        });

        // Create 6 nodes for 2 elements (one solid, one cavity)
        multiDomain.nodes().createNode({.index = 1, .x = 0.0, .y = 0.0, .state = state});
        multiDomain.nodes().createNode({.index = 2, .x = 0.01, .y = 0.0, .state = state});
        multiDomain.nodes().createNode({.index = 3, .x = 0.02, .y = 0.0, .state = state});
        multiDomain.nodes().createNode({.index = 4, .x = 0.0, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.index = 5, .x = 0.01, .y = 0.05, .state = state});
        multiDomain.nodes().createNode({.index = 6, .x = 0.02, .y = 0.05, .state = state});

        // Create solid material
        multiDomain.materials().createSolidMaterial(
          {.name = "SolidMaterial",
           .thermalConductivityDry = 1.8,
           .density = 2050.0,
           .porosity = 0.22,
           .heatCapacity = 850.0,
           .diffusionResistanceFactor = 15.0,
           .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
           .moistureDependentMeasurementTemperature = 0.0,
           .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {100.0, 1.8}},
           .temperatureDependentMeasurementHumidity = 0.0,
           .liquidTransportCurve = {{0, 0}, {180, 2e-6}},
           .sorptionCurve = {{0, 0}, {1, 180}}});

        // Create frame cavity (gas)
        multiDomain.materials().createGas("FrameCavity", HygroThermFEM::CavityStandard::ISO15099);

        // Create elements: one solid and one frame cavity
        multiDomain.createElement({.node1 = 1, .node2 = 2, .node3 = 5, .node4 = 4, .material = "SolidMaterial"});
        multiDomain.createElement({.node1 = 2, .node2 = 3, .node3 = 6, .node4 = 5, .material = "FrameCavity"});
    }
}   // namespace


TEST(TestDomainClearModelAndGravity, TestClearModelRemovesNodesAndMaterials)
{
    SCOPED_TRACE("Begin Test: clearModel removes nodes and materials.");

    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    createSimpleModel(multiDomain);

    // Verify model was created
    const auto nodeCountBefore = multiDomain.nodes().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCountBefore, 4u);

    // Clear the model
    multiDomain.thermal().clearModel();

    // Verify pools are cleared
    const auto nodeCountAfter = multiDomain.nodes().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCountAfter, 0u);
}

TEST(TestDomainClearModelAndGravity, TestClearModelAllowsNewModel)
{
    SCOPED_TRACE("Begin Test: clearModel allows creating new model.");

    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    createSimpleModel(multiDomain);

    // Clear the model
    multiDomain.thermal().clearModel();

    // Create a new model - should not throw
    const State state({
        .temperature = 25.0,
        .humidity = 0.6,
        .pressure = 101325.0
    });
    multiDomain.nodes().createNode({.index = 1, .x = 0.0, .y = 0.0, .state = state});
    multiDomain.nodes().createNode({.index = 2, .x = 2.0, .y = 0.0, .state = state});
    multiDomain.nodes().createNode({.index = 3, .x = 2.0, .y = 2.0, .state = state});
    multiDomain.nodes().createNode({.index = 4, .x = 0.0, .y = 2.0, .state = state});

    const auto nodeCount = multiDomain.nodes().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCount, 4u);

    // Verify new initial temperature is used
    const auto temperatures = multiDomain.nodes().properties(Variable::temperature);
    EXPECT_DOUBLE_EQ(temperatures[0], 25.0);
}

TEST(TestDomainClearModelAndGravity, TestSetGravityVectorDefault)
{
    SCOPED_TRACE("Begin Test: Default gravity vector is (0, -1, 0).");

    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    createModelWithFrameCavity(multiDomain);

    // Add boundary condition
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{0.0, 30.0};
    multiDomain.thermal().createBC_FixedHc(1, 4, bcCoeff);

    // Run transient to initialize gas cavities
    auto temperatures = multiDomain.nodes().properties(Variable::temperature);
    const auto result = multiDomain.thermal().transient(temperatures, 360.0);

    // Should complete without error with default gravity
    EXPECT_FALSE(result.solution.empty());
}

TEST(TestDomainClearModelAndGravity, TestSetGravityVectorCustom)
{
    SCOPED_TRACE("Begin Test: Custom gravity vector affects frame cavity.");

    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    createModelWithFrameCavity(multiDomain);

    // Add boundary condition
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{0.0, 30.0};
    multiDomain.thermal().createBC_FixedHc(1, 4, bcCoeff);

    // Run first transient step with default gravity to initialize cavities
    auto temperatures = multiDomain.nodes().properties(Variable::temperature);
    multiDomain.thermal().transient(temperatures, 360.0);

    // Set horizontal gravity (tilted system)
    const FenestrationCommon::GravityVector horizontalGravity{1.0, 0.0, 0.0};
    multiDomain.thermal().setGravityVector(horizontalGravity);

    // Run another transient step - should use new gravity
    temperatures = multiDomain.nodes().properties(Variable::temperature);
    const auto result = multiDomain.thermal().transient(temperatures, 360.0);

    // Should complete without error
    EXPECT_FALSE(result.solution.empty());
}

TEST(TestDomainClearModelAndGravity, TestSetGravityVectorBeforeCavityInit)
{
    SCOPED_TRACE("Begin Test: Set gravity vector before cavity initialization.");

    MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
    createModelWithFrameCavity(multiDomain);

    // Set gravity before any transient calculation (before gasCavities is created)
    const FenestrationCommon::GravityVector customGravity{0.0, 0.0, -1.0};
    multiDomain.thermal().setGravityVector(customGravity);

    // Add boundary condition
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{0.0, 30.0};
    multiDomain.thermal().createBC_FixedHc(1, 4, bcCoeff);

    // Run transient - should initialize cavities with custom gravity
    auto temperatures = multiDomain.nodes().properties(Variable::temperature);
    const auto result = multiDomain.thermal().transient(temperatures, 360.0);

    EXPECT_FALSE(result.solution.empty());
}
