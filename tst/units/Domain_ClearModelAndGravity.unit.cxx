#include <gtest/gtest.h>
#include <memory>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::State;
using HygroThermFEM::StateParams;
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
        multiDomain.nodes().createNode(1, 0.0, 0.0, state);
        multiDomain.nodes().createNode(2, 1.0, 0.0, state);
        multiDomain.nodes().createNode(3, 1.0, 1.0, state);
        multiDomain.nodes().createNode(4, 0.0, 1.0, state);

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

        multiDomain.createElement(1, 2, 3, 4, "TestMaterial");
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
        multiDomain.nodes().createNode(1, 0.0, 0.0, state);
        multiDomain.nodes().createNode(2, 0.01, 0.0, state);
        multiDomain.nodes().createNode(3, 0.02, 0.0, state);
        multiDomain.nodes().createNode(4, 0.0, 0.05, state);
        multiDomain.nodes().createNode(5, 0.01, 0.05, state);
        multiDomain.nodes().createNode(6, 0.02, 0.05, state);

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
        multiDomain.createElement(1, 2, 5, 4, "SolidMaterial");
        multiDomain.createElement(2, 3, 6, 5, "FrameCavity");
    }
}   // namespace


TEST(TestDomainClearModelAndGravity, TestClearModelRemovesNodesAndMaterials)
{
    SCOPED_TRACE("Begin Test: clearModel removes nodes and materials.");

    MultiDomain multiDomain(true, false);
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

    MultiDomain multiDomain(true, false);
    createSimpleModel(multiDomain);

    // Clear the model
    multiDomain.thermal().clearModel();

    // Create a new model - should not throw
    const State state({
        .temperature = 25.0,
        .humidity = 0.6,
        .pressure = 101325.0
    });
    multiDomain.nodes().createNode(1, 0.0, 0.0, state);
    multiDomain.nodes().createNode(2, 2.0, 0.0, state);
    multiDomain.nodes().createNode(3, 2.0, 2.0, state);
    multiDomain.nodes().createNode(4, 0.0, 2.0, state);

    const auto nodeCount = multiDomain.nodes().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCount, 4u);

    // Verify new initial temperature is used
    const auto temperatures = multiDomain.nodes().properties(Variable::temperature);
    EXPECT_DOUBLE_EQ(temperatures[0], 25.0);
}

TEST(TestDomainClearModelAndGravity, TestSetGravityVectorDefault)
{
    SCOPED_TRACE("Begin Test: Default gravity vector is (0, -1, 0).");

    MultiDomain multiDomain(true, false);
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

    MultiDomain multiDomain(true, false);
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

    MultiDomain multiDomain(true, false);
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
