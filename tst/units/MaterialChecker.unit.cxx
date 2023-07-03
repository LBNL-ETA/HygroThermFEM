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
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
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
    const auto matCheckTransient{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::Transient)};
    EXPECT_EQ(matCheckTransient.size(), 1u);
    EXPECT_EQ(matCheckTransient.isMaterialLibraryCorrect(), false);

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
    const auto matCheckSteadyState{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::SteadyState)};
    EXPECT_EQ(matCheckSteadyState.size(), 1u);
    EXPECT_EQ(matCheckSteadyState.isMaterialLibraryCorrect(), false);
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
    const auto matCheckTransient{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::Transient)};
    EXPECT_EQ(matCheckTransient.size(), 1u);
    EXPECT_EQ(matCheckTransient.isMaterialLibraryCorrect(), false);

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
    const auto matCheckSteadyState{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::SteadyState)};
    EXPECT_EQ(matCheckSteadyState.size(), 1u);
    EXPECT_EQ(matCheckSteadyState.isMaterialLibraryCorrect(), false);
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

TEST_F(TestMaterialChecker, TestExample_3)
{
    SCOPED_TRACE("Test for material missing properties in case of only moisture simulation.");

    // Note that material is incomplete
    MaterialPool::Instance().createSolidMaterial("Test material");

    const bool SimulateThermal{false};
    const bool SimulateMoisture{true};
    HygroThermFEM::MultiDomain domain{SimulateThermal, SimulateMoisture};

    // Check for transient simulation
    const auto matCheckTransient{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::Transient)};
    EXPECT_EQ(matCheckTransient.size(), 1u);
    EXPECT_EQ(matCheckTransient.isMaterialLibraryCorrect(), false);

    // Series of test means following. Material property is missing for given simulation properties
    // and therefore it is set to true. If material property is not needed for given simulation
    // properties then it is set to false.
    EXPECT_EQ(matCheckTransient[0].Density, false);
    EXPECT_EQ(matCheckTransient[0].Emissivity, false);
    EXPECT_EQ(matCheckTransient[0].Porosity, true);
    EXPECT_EQ(matCheckTransient[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityDry, false);
    EXPECT_EQ(matCheckTransient[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckTransient[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    // Check for steady-state simulation
    const auto matCheckSteadyState{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::SteadyState)};
    EXPECT_EQ(matCheckSteadyState.size(), 1u);
    EXPECT_EQ(matCheckSteadyState.isMaterialLibraryCorrect(), false);
    EXPECT_EQ(matCheckSteadyState[0].Density, false);
    EXPECT_EQ(matCheckSteadyState[0].Emissivity, false);
    EXPECT_EQ(matCheckSteadyState[0].Porosity, true);
    EXPECT_EQ(matCheckSteadyState[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckSteadyState[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    MaterialPool::Instance().clear();
}

TEST_F(TestMaterialChecker, TestExample_4)
{
    SCOPED_TRACE("Test for material missing properties with variable simulation parameters "
                 "(exclude water liquid transportation).");

    // Note that material is incomplete
    MaterialPool::Instance().createSolidMaterial("Test material");

    const bool SimulateThermal{true};
    const bool SimulateMoisture{true};
    HygroThermFEM::MultiDomain domain{SimulateThermal, SimulateMoisture};

    // Simulation properties will have turned off liquid transportation and capillary conduction in
    // order to test need for liquid transportation curve.
    const auto excludeWaterLiquidTransportation{true};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{true};
    const auto excludeVaporDiffusionConduction{false};
    const auto thermalConductivityMoistureAndTemperatureDependent{false};

    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    // Check for transient simulation
    const auto matCheckTransient{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::Transient)};
    EXPECT_EQ(matCheckTransient.size(), 1u);
    EXPECT_EQ(matCheckTransient.isMaterialLibraryCorrect(), false);

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
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationSuction, false);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    // Check for steady-state simulation
    const auto matCheckSteadyState{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::SteadyState)};
    EXPECT_EQ(matCheckSteadyState.size(), 1u);
    EXPECT_EQ(matCheckSteadyState.isMaterialLibraryCorrect(), false);
    EXPECT_EQ(matCheckSteadyState[0].Density, false);
    EXPECT_EQ(matCheckSteadyState[0].Emissivity, true);
    EXPECT_EQ(matCheckSteadyState[0].Porosity, true);
    EXPECT_EQ(matCheckSteadyState[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckSteadyState[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckSteadyState[0].MoistureStorageFunction, false);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationSuction, false);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    MaterialPool::Instance().clear();
}

TEST_F(TestMaterialChecker, TestExample_5)
{
    SCOPED_TRACE("Test for material missing properties with variable simulation parameters "
                 "(water vapor diffusion resistance factor not needed).");

    // Note that material is incomplete
    MaterialPool::Instance().createSolidMaterial("Test material");

    const bool SimulateThermal{true};
    const bool SimulateMoisture{false};
    HygroThermFEM::MultiDomain domain{SimulateThermal, SimulateMoisture};

    // Simulation properties will have turned off liquid transportation and capillary conduction in
    // order to test need for liquid transportation curve.
    const auto excludeWaterLiquidTransportation{false};
    const auto excludeHeatOfEvaporation{true};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{false};
    const auto thermalConductivityMoistureAndTemperatureDependent{false};

    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    // Check for transient simulation
    const auto matCheckTransient{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::Transient)};
    EXPECT_EQ(matCheckTransient.size(), 1u);
    EXPECT_EQ(matCheckTransient.isMaterialLibraryCorrect(), false);

    // Series of test means following. Material property is missing for given simulation properties
    // and therefore it is set to true. If material property is not needed for given simulation
    // properties then it is set to false.
    EXPECT_EQ(matCheckTransient[0].Density, true);
    EXPECT_EQ(matCheckTransient[0].Emissivity, true);
    EXPECT_EQ(matCheckTransient[0].Porosity, false);
    EXPECT_EQ(matCheckTransient[0].SpecificHeatCapacityDry, true);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckTransient[0].WaterVaporDiffusionResistanceFactor, false);
    EXPECT_EQ(matCheckTransient[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckTransient[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckTransient[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    // Check for steady-state simulation
    const auto matCheckSteadyState{
      checkForMaterialsValidity(domain, HygroThermFEM::SimulationType::SteadyState)};
    EXPECT_EQ(matCheckSteadyState.size(), 1u);
    EXPECT_EQ(matCheckSteadyState.isMaterialLibraryCorrect(), false);
    EXPECT_EQ(matCheckSteadyState[0].Density, false);
    EXPECT_EQ(matCheckSteadyState[0].Emissivity, true);
    EXPECT_EQ(matCheckSteadyState[0].Porosity, false);
    EXPECT_EQ(matCheckSteadyState[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckSteadyState[0].WaterVaporDiffusionResistanceFactor, false);
    EXPECT_EQ(matCheckSteadyState[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationSuction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityMoistureAndTemperatureDependent, false);

    MaterialPool::Instance().clear();
}