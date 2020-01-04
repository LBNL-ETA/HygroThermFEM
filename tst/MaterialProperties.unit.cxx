#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

class TestMaterialProperties : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        HygroThermFEM::NodePool::Instance().clear();
        HygroThermFEM::MaterialPool::Instance().clear();
    }
};

TEST_F(TestMaterialProperties, PropertiesExistence)
{
    SCOPED_TRACE("Material properties.");

    // Material Properties (Cottaer Sandstone)
    const double thermalConductivityDry{1.8};
    const double density{2050.0};
    const double porosity{0.22};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    const double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                              {27, 1E-8},
                                                                              {45, 1.1E-8},
                                                                              {90, 2E-8},
                                                                              {126, 3.5E-8},
                                                                              {144, 5E-8},
                                                                              {162, 1E-7},
                                                                              {171, 2E-7},
                                                                              {180, 7E-7}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                            {0.5, 5.3},
                                                                            {0.65, 8.4},
                                                                            {0.8, 12},
                                                                            {0.93, 17},
                                                                            {0.95, 25},
                                                                            {0.99, 63},
                                                                            {0.995, 83},
                                                                            {0.999, 120},
                                                                            {1, 180}};

    auto & material =
      HygroThermFEM::MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone");

    EXPECT_EQ(false, material.hasThermalConductivityDry());
    EXPECT_EQ(false, material.hasPorosity());
    EXPECT_EQ(false, material.hasDensity());
    EXPECT_EQ(false, material.hasHeatCapacity());
    EXPECT_EQ(false, material.hasDiffusionResistanceFactor());
    EXPECT_EQ(false, material.hasLiquidTransportationCurve());
    EXPECT_EQ(false, material.hasSorptionCurve());
    EXPECT_EQ(false, material.hasThermalConductivityMoistureAndTemperatureDependent());

    material.setThermalConductivity(thermalConductivityDry);
    material.setPorosity(porosity);
    material.setDensity(density);
    material.setHeatCapacity(specificHeatCapacityDry);
    material.setDiffusionResistanceFactor(diffusionResistanceFactor);
    material.setLiquidTransportationCurve(liquidTransportationCurve);
    material.setSorptionCurve(moistureStorageFunction);
    material.setThermalConductivityMoistureAndTemperatureDependent(
      thermalConductivityMoistureDependent,
      thermalConductivityMeasuredAtTemperature,
      thermalConductivityTemperatureDependent,
      thermalConductivityMeasuredAtHumidity);

    EXPECT_EQ(true, material.hasThermalConductivityDry());
    EXPECT_EQ(true, material.hasPorosity());
    EXPECT_EQ(true, material.hasDensity());
    EXPECT_EQ(true, material.hasHeatCapacity());
    EXPECT_EQ(true, material.hasDiffusionResistanceFactor());
    EXPECT_EQ(true, material.hasLiquidTransportationCurve());
    EXPECT_EQ(true, material.hasSorptionCurve());
    EXPECT_EQ(true, material.hasThermalConductivityMoistureAndTemperatureDependent());
}
