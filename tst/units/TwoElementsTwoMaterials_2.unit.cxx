#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::WaterContent;
using HygroThermFEM::MultiDomain;

/// This is example of water content calculation in nodes that are shared between two elements
/// with different material. In this case, one of the elements is triangular which will lead
/// to different water content solution at the shared nodes.
TEST(TwoElementsTwoMaterials_2, NodeInTwoMaterials)
{
    SCOPED_TRACE("Begin Test: Node as part of two elements that have different material and one "
                 "element is triangular.");

    constexpr auto temperature = 10.0;
    constexpr auto humidity = 0.8;
    constexpr auto pressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::State state({
        .temperature = temperature,
        .humidity = humidity,
        .pressure = pressure,
        .liquidPercent = liquidPercent
    });

    MultiDomain multiDomain(false, false);

    multiDomain.nodes().createNode(1, 0, 1, state);
    multiDomain.nodes().createNode(2, 1, 0, state);
    multiDomain.nodes().createNode(3, 1, 1, state);
    multiDomain.nodes().createNode(4, 2, 0, state);
    multiDomain.nodes().createNode(5, 2, 1, state);

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

    multiDomain.createElement(1, 1, 2, 3, material1.name());
    multiDomain.createElement(2, 4, 5, 3, material2.name());

    auto iceContent = multiDomain.property(HygroThermFEM::Variable::ice);
    auto vaporContent = multiDomain.property(HygroThermFEM::Variable::vapor);
    auto liquidContent = multiDomain.property(HygroThermFEM::Variable::liquid);

    /// Test water content in node number 2 (material 2 will have more influence)
    EXPECT_NEAR(iceContent[1], 0, 1e-6);
    EXPECT_NEAR(vaporContent[1], 0.000903, 1e-6);
    EXPECT_NEAR(liquidContent[1], 60.66576371, 1e-6);

    /// Node number 3 should have different water content from node number 2
    /// Influence of materials in this node is identical.
    EXPECT_NEAR(iceContent[2], 0, 1e-6);
    EXPECT_NEAR(vaporContent[2], 0.001062, 1e-6);
    EXPECT_NEAR(liquidContent[2], 48.498938, 1e-6);
}
