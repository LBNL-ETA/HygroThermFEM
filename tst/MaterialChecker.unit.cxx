#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class TestMaterialChecker : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(TestMaterialChecker, TestExample_1)
{
    SCOPED_TRACE(
      "Test for material missing properties in case of thermal and moisture simulations are on.");

    // Note that material is incomplete
    MaterialPool::Instance().createSolidMaterial("Test material");

    const bool SimulateThermal{true};
    const bool SimulateMoisture{true};
    HygroThermFEM::MultiDomain domain{SimulateThermal, SimulateMoisture};

    // Check for transient simulation
    const auto matCheckTransient = domain.checkMaterialsForTransientSimulation();
    EXPECT_EQ(matCheckTransient.size(), 1);

    // Series of test means following. Material property is missing for given simulation properties
    // and therefore it is set to true. If material property is not needed for given simulation
    // properties then it is set to false.
    EXPECT_EQ(matCheckTransient[0].Density, true);
    EXPECT_EQ(matCheckTransient[0].Emissivity, true);
    EXPECT_EQ(matCheckTransient[0].Porosity, true);
    EXPECT_EQ(matCheckTransient[0].SpecificHeatCapacityDry, true);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckTransient[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckTransient[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    // Check for steady-state simulation
    const auto matCheckSteadyState = domain.checkMaterialsForSteadyStateSimulation();
    EXPECT_EQ(matCheckSteadyState.size(), 1);
    EXPECT_EQ(matCheckSteadyState[0].Density, false);
    EXPECT_EQ(matCheckSteadyState[0].Emissivity, true);
    EXPECT_EQ(matCheckSteadyState[0].Porosity, true);
    EXPECT_EQ(matCheckSteadyState[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckSteadyState[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckSteadyState[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    MaterialPool::Instance().clear();
}

TEST_F(TestMaterialChecker, TestExample_2)
{
    SCOPED_TRACE("Test for material missing properties in case of only thermal simulation.");

    // Note that material is incomplete
    MaterialPool::Instance().createSolidMaterial("Test material");

    const bool SimulateThermal{true};
    const bool SimulateMoisture{false};
    HygroThermFEM::MultiDomain domain{SimulateThermal, SimulateMoisture};

    // Check for transient simulation
    const auto matCheckTransient = domain.checkMaterialsForTransientSimulation();
    EXPECT_EQ(matCheckTransient.size(), 1);

    // Series of test means following. Material property is missing for given simulation properties
    // and therefore it is set to true. If material property is not needed for given simulation
    // properties then it is set to false.
    EXPECT_EQ(matCheckTransient[0].Density, true);
    EXPECT_EQ(matCheckTransient[0].Emissivity, true);
    EXPECT_EQ(matCheckTransient[0].Porosity, false);
    EXPECT_EQ(matCheckTransient[0].SpecificHeatCapacityDry, true);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckTransient[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckTransient[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    // Check for steady-state simulation
    const auto matCheckSteadyState = domain.checkMaterialsForSteadyStateSimulation();
    EXPECT_EQ(matCheckSteadyState.size(), 1);
    EXPECT_EQ(matCheckSteadyState[0].Density, false);
    EXPECT_EQ(matCheckSteadyState[0].Emissivity, true);
    EXPECT_EQ(matCheckSteadyState[0].Porosity, false);
    EXPECT_EQ(matCheckSteadyState[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckSteadyState[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckSteadyState[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    MaterialPool::Instance().clear();
}
