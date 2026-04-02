#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestMaterialProperties, PropertiesExistence)
{
    SCOPED_TRACE("Material properties.");

    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
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

    HygroThermFEM::MultiDomain multiDomain;

    auto & material = multiDomain.materials().createSolidMaterial("Cottaer Sandstone");

    EXPECT_EQ(false, material.hasThermalConductivityDry());
    EXPECT_EQ(false, material.hasPorosity());
    EXPECT_EQ(false, material.hasDensity());
    EXPECT_EQ(false, material.hasHeatCapacity());
    EXPECT_EQ(false, material.hasDiffusionResistanceFactor());
    EXPECT_EQ(false, material.hasLiquidTransportationCurve());
    EXPECT_EQ(false, material.hasSorptionCurve());
    EXPECT_EQ(false, material.hasThermalConductivityMoistureAndTemperatureDependent());

    material.setThermalConductivity(1.8);
    material.setPorosity(0.22);
    material.setDensity(2050.0);
    material.setHeatCapacity(850.0);
    material.setDiffusionResistanceFactor(15.0);
    material.setLiquidTransportationCurve(liquidTransportationCurve);
    material.setSorptionCurve(moistureStorageFunction);
    material.setThermalConductivityMoistureAndTemperatureDependent(
      thermalConductivityMoistureDependent,
      0,
      thermalConductivityTemperatureDependent,
      0);

    EXPECT_EQ(true, material.hasThermalConductivityDry());
    EXPECT_EQ(true, material.hasPorosity());
    EXPECT_EQ(true, material.hasDensity());
    EXPECT_EQ(true, material.hasHeatCapacity());
    EXPECT_EQ(true, material.hasDiffusionResistanceFactor());
    EXPECT_EQ(true, material.hasLiquidTransportationCurve());
    EXPECT_EQ(true, material.hasSorptionCurve());
    EXPECT_EQ(true, material.hasThermalConductivityMoistureAndTemperatureDependent());
}
