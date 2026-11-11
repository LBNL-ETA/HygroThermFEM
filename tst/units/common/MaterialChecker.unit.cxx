#include <memory>
#include <vector>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::State;

namespace
{
    //! Physics in which nothing but latent heat of fusion reads a material's water:
    //! capillary conduction, heat of evaporation and vapour diffusion are all off.
    HygroThermFEM::PhysicsOptions onlyFusionReadsWater(const bool modelFusion)
    {
        return {.excludeHeatOfEvaporation = true,
                .excludeCapillaryConduction = true,
                .excludeVaporDiffusionConduction = true,
                .excludeLatentHeatOfFusion = !modelFusion};
    }

    const HygroThermFEM::MultiDomainParams heatOnly{.performThermal = true, .performMoisture = false};

    const std::vector<HygroThermFEM::SimulationType> bothModes{HygroThermFEM::SimulationType::Transient,
                                                              HygroThermFEM::SimulationType::SteadyState};
}   // namespace

TEST(MaterialCheckerFusion, FusionNeedsSorptionCurveAndPorosity)
{
    HygroThermFEM::MultiDomain multiDomain(heatOnly);
    multiDomain.setPhysicsOptions(onlyFusionReadsWater(true));
    multiDomain.materials().createSolidMaterial("No water data");

    for(const auto mode : bothModes)
    {
        const auto check{multiDomain.checkForMaterialsValidity(mode)};
        ASSERT_EQ(1u, check.size());
        EXPECT_TRUE(check[0].MoistureStorageFunction);
        EXPECT_TRUE(check[0].Porosity);
    }
}

TEST(MaterialCheckerFusion, WithoutFusionNeitherIsNeeded)
{
    HygroThermFEM::MultiDomain multiDomain(heatOnly);
    multiDomain.setPhysicsOptions(onlyFusionReadsWater(false));
    multiDomain.materials().createSolidMaterial("No water data");

    for(const auto mode : bothModes)
    {
        const auto check{multiDomain.checkForMaterialsValidity(mode)};
        ASSERT_EQ(1u, check.size());
        EXPECT_FALSE(check[0].MoistureStorageFunction);
        EXPECT_FALSE(check[0].Porosity);
    }
}

TEST(MaterialCheckerFusion, CurveWithoutPorosityIsReported)
{
    HygroThermFEM::MultiDomain multiDomain(heatOnly);
    multiDomain.setPhysicsOptions(onlyFusionReadsWater(true));
    auto & material{multiDomain.materials().createSolidMaterial("Curve only")};
    material.setSorptionCurve({{0.0, 0.0}, {0.5, 10.0}, {1.0, 100.0}});

    const auto check{multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient)};
    ASSERT_EQ(1u, check.size());
    EXPECT_FALSE(check[0].MoistureStorageFunction);
    EXPECT_TRUE(check[0].Porosity);
}

//! A material that holds no water says so with a zero curve and zero porosity, and that
//! satisfies the check: the data is there, and it is zero.
TEST(MaterialCheckerFusion, ZeroCurveAndZeroPorosityAreData)
{
    HygroThermFEM::MultiDomain multiDomain(heatOnly);
    multiDomain.setPhysicsOptions(onlyFusionReadsWater(true));
    auto & material{multiDomain.materials().createSolidMaterial("Aluminum")};
    material.setSorptionCurve({{0.0, 0.0}, {1.0, 0.0}});
    material.setPorosity(0.0);

    const auto check{multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient)};
    ASSERT_EQ(1u, check.size());
    EXPECT_FALSE(check[0].MoistureStorageFunction);
    EXPECT_FALSE(check[0].Porosity);
}

class TestMaterialChecker : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        HygroThermFEM::SimulationProperties::Instance().resetCalculationParameters();
    }
};

TEST_F(TestMaterialChecker, TestExample_1)
{
    SCOPED_TRACE(
      "Test for material missing properties in case of thermal and moisture simulations are on.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

    // Note that material is incomplete
    multiDomain.materials().createSolidMaterial("Test material");

    // Check for transient simulation
    const auto matCheckTransient = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient);
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
    const auto matCheckSteadyState = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::SteadyState);
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
}

TEST_F(TestMaterialChecker, TestExample_2)
{
    SCOPED_TRACE("Test for material missing properties in case of only thermal simulation.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Note that material is incomplete
    multiDomain.materials().createSolidMaterial("Test material");

    // Check for transient simulation
    const auto matCheckTransient = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient);
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
    const auto matCheckSteadyState = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::SteadyState);
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
}

TEST_F(TestMaterialChecker, TestExample_3)
{
    SCOPED_TRACE("Test for material missing properties in case of only moisture simulation.");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    // Note that material is incomplete
    multiDomain.materials().createSolidMaterial("Test material");

    // Check for transient simulation
    const auto matCheckTransient = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient);
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
    const auto matCheckSteadyState = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::SteadyState);
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
}

TEST_F(TestMaterialChecker, TestExample_4)
{
    SCOPED_TRACE("Test for material missing properties with variable simulation parameters "
                 "(exclude water liquid transportation).");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

    // Note that material is incomplete
    multiDomain.materials().createSolidMaterial("Test material");

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
    const auto matCheckTransient = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient);
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
    const auto matCheckSteadyState = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::SteadyState);
    EXPECT_EQ(matCheckSteadyState.size(), 1u);
    EXPECT_EQ(matCheckSteadyState.isMaterialLibraryCorrect(), false);
    EXPECT_EQ(matCheckSteadyState[0].Density, false);
    EXPECT_EQ(matCheckSteadyState[0].Emissivity, true);
    EXPECT_EQ(matCheckSteadyState[0].Porosity, true);
    EXPECT_EQ(matCheckSteadyState[0].SpecificHeatCapacityDry, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityDry, true);
    EXPECT_EQ(matCheckSteadyState[0].WaterVaporDiffusionResistanceFactor, true);
    EXPECT_EQ(matCheckSteadyState[0].MoistureStorageFunction, true);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationSuction, false);
    EXPECT_EQ(matCheckSteadyState[0].LiquidTransportationRedistribution, false);
    EXPECT_EQ(matCheckSteadyState[0].ThermalConductivityMoistureAndTemperatureDependent, false);
}

TEST_F(TestMaterialChecker, TestExample_5)
{
    SCOPED_TRACE("Test for material missing properties with variable simulation parameters "
                 "(water vapor diffusion resistance factor not needed).");

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Note that material is incomplete
    multiDomain.materials().createSolidMaterial("Test material");

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
    const auto matCheckTransient = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::Transient);
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
    const auto matCheckSteadyState = multiDomain.checkForMaterialsValidity(HygroThermFEM::SimulationType::SteadyState);
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
}