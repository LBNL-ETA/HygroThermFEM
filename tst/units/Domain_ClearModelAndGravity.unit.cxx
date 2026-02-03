#include <gtest/gtest.h>
#include <memory>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;
using HygroThermFEM::ThermalDomain;
using HygroThermFEM::MultiDomain;
using HygroThermFEM::Variable;

class TestDomainClearModelAndGravity : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }

    // Helper to create a simple 2x2 model with one element
    void createSimpleModel(MultiDomain & multiDomain)
    {
        const State state(20.0, 0.5, 101325.0);

        // Create 4 nodes for a single quad element
        NodePool::Instance().createNode(1, 0.0, 0.0, state);
        NodePool::Instance().createNode(2, 1.0, 0.0, state);
        NodePool::Instance().createNode(3, 1.0, 1.0, state);
        NodePool::Instance().createNode(4, 0.0, 1.0, state);

        // Create material
        constexpr double thermalConductivityDry{1.0};
        constexpr double density{1000.0};
        constexpr double porosity{0.1};
        constexpr double specificHeatCapacityDry{1000.0};
        constexpr double diffusionResistanceFactor{10.0};
        const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
          {0.0, 1.0}, {180, 1.0}};
        constexpr double thermalConductivityMeasuredAtTemperature{0};
        const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
          {0.0, 1.0}, {100.0, 1.0}};
        constexpr double thermalConductivityMeasuredAtHumidity{0};
        const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0}, {180, 1e-6}};
        const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0}, {1, 180}};

        MaterialPool::Instance().createSolidMaterial("TestMaterial",
                                                     thermalConductivityDry,
                                                     density,
                                                     porosity,
                                                     specificHeatCapacityDry,
                                                     diffusionResistanceFactor,
                                                     thermalConductivityMoistureDependent,
                                                     thermalConductivityMeasuredAtTemperature,
                                                     thermalConductivityTemperatureDependent,
                                                     thermalConductivityMeasuredAtHumidity,
                                                     liquidTransportationCurve,
                                                     moistureStorageFunction);

        multiDomain.createElement(1, 2, 3, 4, "TestMaterial");
    }

    // Helper to create model with frame cavity for gravity tests
    void createModelWithFrameCavity(MultiDomain & multiDomain)
    {
        const State state(20.0, 0.0, 101325.0);

        // Create 6 nodes for 2 elements (one solid, one cavity)
        NodePool::Instance().createNode(1, 0.0, 0.0, state);
        NodePool::Instance().createNode(2, 0.01, 0.0, state);
        NodePool::Instance().createNode(3, 0.02, 0.0, state);
        NodePool::Instance().createNode(4, 0.0, 0.05, state);
        NodePool::Instance().createNode(5, 0.01, 0.05, state);
        NodePool::Instance().createNode(6, 0.02, 0.05, state);

        // Create solid material
        constexpr double thermalConductivityDry{1.8};
        constexpr double density{2050.0};
        constexpr double porosity{0.22};
        constexpr double specificHeatCapacityDry{850.0};
        constexpr double diffusionResistanceFactor{15.0};
        const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
          {0.0, 1.8}, {180, 1.8}};
        const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
          {0.0, 1.8}, {100.0, 1.8}};
        const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0}, {180, 2e-6}};
        const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0}, {1, 180}};

        MaterialPool::Instance().createSolidMaterial("SolidMaterial",
                                                     thermalConductivityDry,
                                                     density,
                                                     porosity,
                                                     specificHeatCapacityDry,
                                                     diffusionResistanceFactor,
                                                     thermalConductivityMoistureDependent,
                                                     0.0,
                                                     thermalConductivityTemperatureDependent,
                                                     0.0,
                                                     liquidTransportationCurve,
                                                     moistureStorageFunction);

        // Create frame cavity (gas)
        MaterialPool::Instance().createGas("FrameCavity", HygroThermFEM::CavityStandard::ISO15099);

        // Create elements: one solid and one frame cavity
        multiDomain.createElement(1, 2, 5, 4, "SolidMaterial");
        multiDomain.createElement(2, 3, 6, 5, "FrameCavity");
    }
};

TEST_F(TestDomainClearModelAndGravity, TestClearModelRemovesNodesAndMaterials)
{
    SCOPED_TRACE("Begin Test: clearModel removes nodes and materials.");

    MultiDomain multiDomain(true, false);
    auto & domain = multiDomain.thermal();
    createSimpleModel(multiDomain);

    // Verify model was created
    const auto nodeCountBefore = NodePool::Instance().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCountBefore, 4u);

    // Clear the model
    domain.clearModel();

    // Verify pools are cleared
    const auto nodeCountAfter = NodePool::Instance().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCountAfter, 0u);
}

TEST_F(TestDomainClearModelAndGravity, TestClearModelAllowsNewModel)
{
    SCOPED_TRACE("Begin Test: clearModel allows creating new model.");

    MultiDomain multiDomain(true, false);
    auto & domain = multiDomain.thermal();
    createSimpleModel(multiDomain);

    // Clear the model
    domain.clearModel();

    // Create a new model - should not throw
    const State state(25.0, 0.6, 101325.0);
    NodePool::Instance().createNode(1, 0.0, 0.0, state);
    NodePool::Instance().createNode(2, 2.0, 0.0, state);
    NodePool::Instance().createNode(3, 2.0, 2.0, state);
    NodePool::Instance().createNode(4, 0.0, 2.0, state);

    const auto nodeCount = NodePool::Instance().properties(Variable::temperature).size();
    EXPECT_EQ(nodeCount, 4u);

    // Verify new initial temperature is used
    const auto temperatures = NodePool::Instance().properties(Variable::temperature);
    EXPECT_DOUBLE_EQ(temperatures[0], 25.0);
}

TEST_F(TestDomainClearModelAndGravity, TestSetGravityVectorDefault)
{
    SCOPED_TRACE("Begin Test: Default gravity vector is (0, -1, 0).");

    MultiDomain multiDomain(true, false);
    auto & domain = multiDomain.thermal();
    createModelWithFrameCavity(multiDomain);

    // Add boundary condition
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{0.0, 30.0};
    domain.createBC_FixedHc(1, 4, bcCoeff);

    // Run transient to initialize gas cavities
    auto temperatures = NodePool::Instance().properties(Variable::temperature);
    const auto result = domain.transient(temperatures, 360.0);

    // Should complete without error with default gravity
    EXPECT_FALSE(result.solution.empty());
}

TEST_F(TestDomainClearModelAndGravity, TestSetGravityVectorCustom)
{
    SCOPED_TRACE("Begin Test: Custom gravity vector affects frame cavity.");

    MultiDomain multiDomain(true, false);
    auto & domain = multiDomain.thermal();
    createModelWithFrameCavity(multiDomain);

    // Add boundary condition
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{0.0, 30.0};
    domain.createBC_FixedHc(1, 4, bcCoeff);

    // Run first transient step with default gravity to initialize cavities
    auto temperatures = NodePool::Instance().properties(Variable::temperature);
    domain.transient(temperatures, 360.0);

    // Set horizontal gravity (tilted system)
    const FenestrationCommon::GravityVector horizontalGravity{1.0, 0.0, 0.0};
    domain.setGravityVector(horizontalGravity);

    // Run another transient step - should use new gravity
    temperatures = NodePool::Instance().properties(Variable::temperature);
    const auto result = domain.transient(temperatures, 360.0);

    // Should complete without error
    EXPECT_FALSE(result.solution.empty());
}

TEST_F(TestDomainClearModelAndGravity, TestSetGravityVectorBeforeCavityInit)
{
    SCOPED_TRACE("Begin Test: Set gravity vector before cavity initialization.");

    MultiDomain multiDomain(true, false);
    auto & domain = multiDomain.thermal();
    createModelWithFrameCavity(multiDomain);

    // Set gravity before any transient calculation (before gasCavities is created)
    const FenestrationCommon::GravityVector customGravity{0.0, 0.0, -1.0};
    domain.setGravityVector(customGravity);

    // Add boundary condition
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{0.0, 30.0};
    domain.createBC_FixedHc(1, 4, bcCoeff);

    // Run transient - should initialize cavities with custom gravity
    auto temperatures = NodePool::Instance().properties(Variable::temperature);
    const auto result = domain.transient(temperatures, 360.0);

    EXPECT_FALSE(result.solution.empty());
}
