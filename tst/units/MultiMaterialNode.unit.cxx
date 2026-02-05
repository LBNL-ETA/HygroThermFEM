#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::WaterContent;
using HygroThermFEM::MultiDomain;

TEST(MultiMaterialNode, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Node with multiple materials in it.");

    constexpr auto temperature = 10;
    constexpr auto humidity = 0.8;
    constexpr auto pressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::State state({
        .temperature = temperature,
        .humidity = humidity,
        .pressure = pressure,
        .liquidPercent = liquidPercent
    });

    MultiDomain multiDomain;

    auto node1 = multiDomain.nodes().createNode(1, 0, 0, state);

    // Material Properties (Cottaer Sandstone)
    const auto & material1 = multiDomain.materials().createSolidMaterial({
        .name = "Cottaer Sandstone",
        .thermalConductivityDry = 1.8,
        .density = 2050.0,
        .porosity = 0.22,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 1.8}},
        .moistureDependentMeasurementTemperature = 0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 1.8}},
        .temperatureDependentMeasurementHumidity = 0,
        .liquidTransportCurve = {{0, 0},
                                 {27, 1E-8},
                                 {45, 1.1E-8},
                                 {90, 2E-8},
                                 {126, 3.5E-8},
                                 {144, 5E-8},
                                 {162, 1E-7},
                                 {171, 2E-7},
                                 {180, 7E-7}},
        .sorptionCurve = {{0, 0},
                          {0.5, 5.3},
                          {0.65, 8.4},
                          {0.8, 12},
                          {0.93, 17},
                          {0.95, 25},
                          {0.99, 63},
                          {0.995, 83},
                          {0.999, 120},
                          {1, 180}}
    });

    // Material Properties (Concrete, w/c=0.5)
    const auto & material2 = multiDomain.materials().createSolidMaterial({
        .name = "Concrete, w/c=0.5",
        .thermalConductivityDry = 1.6,
        .density = 2300,
        .porosity = 0.18,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 92,
        .thermalConductivityMoistureDependent = {{0.0, 1.6}, {150, 1.6}},
        .moistureDependentMeasurementTemperature = 0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.6}, {1, 1.6}},
        .temperatureDependentMeasurementHumidity = 0,
        .liquidTransportCurve = {{0, 0},
                                 {72, 7.4E-11},
                                 {85, 2.5E-10},
                                 {100, 1E-9},
                                 {118, 1.2E-9},
                                 {150, 1.2e-9}},
        .sorptionCurve = {{0, 0},
                          {0.05, 27},
                          {0.1, 32},
                          {0.15, 34},
                          {0.2, 35},
                          {0.3, 37},
                          {0.4, 40},
                          {0.5, 48},
                          {0.6, 58},
                          {0.7, 72},
                          {0.8, 85},
                          {0.9, 100},
                          {0.95, 118},
                          {1, 150}}
    });

    node1.assignMaterial(material1, 0.5);
    node1.assignMaterial(material2, 0.5);

    EXPECT_NEAR(node1.property(HygroThermFEM::Variable::ice), 0, 1e-6);
    EXPECT_NEAR(node1.property(HygroThermFEM::Variable::vapor), 0.001062, 1e-6);
    EXPECT_NEAR(node1.property(HygroThermFEM::Variable::liquid), 48.498938, 1e-6);
}
