#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "SlabCreator.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ThermalConductivityDependent_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        constexpr auto relaxationParameter{0.8};
        constexpr auto errorTolerance{1e-5};
        constexpr auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ThermalConductivityDependent_1, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{true};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain;

    auto params = TestHelper::CottaerSandstone();
    params.thermalConductivityMoistureDependent = {{0.0, 1.8}, {180, 2.5}};
    params.thermalConductivityTemperatureDependent = {{0.0, 1.8}, {1, 3.1}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    constexpr HygroThermFEM::State state({
        .temperature = 0.0,
        .humidity = 0.99,
        .pressure = 101325.0,
        .liquidPercent = 1.0
    });

    TestHelper::SlabBuilder(multiDomain)
        .gridXCoordinates({0, 0.05, 0.1})
        .height(0.05)
        .material(material.name())
        .state(state)
        .startCorner(TestHelper::StartCorner::BottomRight)
        .direction(TestHelper::Direction::CounterClockwise)
        .build();

    /// Create Boundary Conditions
    constexpr auto hc = 5.0;
    constexpr auto airTemperature = 10.0;
    constexpr auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff);
    multiDomain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

    const auto results = multiDomain.transientMultiStep(dTime, nSteps);

    const std::vector<std::vector<double>> correctWaterContentSolution{{60.6377647, 60.6377647, 62.7016201, 62.7016201, 60.6377647, 60.6377647},
 {58.4639576, 58.4639576, 62.1994746, 62.1994746, 58.4639576, 58.4639576},
 {56.4201644, 56.4201644, 61.5590086, 61.5590086, 56.4201644, 56.4201644},
 {54.4628438, 54.4628438, 60.8262478, 60.8262478, 54.4628438, 54.4628438},
 {52.5606379, 52.5606379, 60.0342318, 60.0342318, 52.5606379, 52.5606379},
 {50.690633, 50.690633, 59.2072722, 59.2072722, 50.690633, 50.690633},
 {48.8355903, 48.8355903, 58.3637879, 58.3637879, 48.8355903, 48.8355903},
 {46.9820692, 46.9820692, 57.5182081, 57.5182081, 46.9820692, 46.9820692},
 {45.1191431, 45.1191431, 56.682286, 56.682286, 45.1191431, 45.1191431},
 {43.2563263, 43.2563263, 55.8471432, 55.8471432, 43.2563263, 43.2563263},
 {41.4115536, 41.4115536, 54.9953012, 54.9953012, 41.4115536, 41.4115536},
 {39.5775178, 39.5775178, 54.1344484, 54.1344484, 39.5775178, 39.5775178},
 {37.7481759, 37.7481759, 53.2709437, 53.2709437, 37.7481759, 37.7481759},
 {35.9184663, 35.9184663, 52.4101124, 52.4101124, 35.9184663, 35.9184663},
 {34.084102, 34.084102, 51.5564629, 51.5564629, 34.084102, 34.084102},
 {32.2413986, 32.2413986, 50.7138673, 50.7138673, 32.2413986, 32.2413986},
 {30.3871429, 30.3871429, 49.8856996, 49.8856996, 30.3871429, 30.3871429},
 {28.5184915, 28.5184915, 49.0749432, 49.0749432, 28.5184915, 28.5184915},
 {26.5948489, 26.5948489, 48.3224296, 48.3224296, 26.5948489, 26.5948489},
 {24.5750913, 24.5750913, 47.6713394, 47.6713394, 24.5750913, 24.5750913},
 {22.5851611, 22.5851611, 47.0018361, 47.0018361, 22.5851611, 22.5851611},
 {20.6290346, 20.6290346, 46.3089648, 46.3089648, 20.6290346, 20.6290346},
 {18.6967376, 18.6967376, 45.6020798, 45.6020798, 18.6967376, 18.6967376},
 {16.8419769, 16.8419769, 44.8389013, 44.8389013, 16.8419769, 16.8419769}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<std::vector<double>> correctTemperatureSolution{{0.117003514, 0.117003514, 0.0805798249, 0.0805798249, 0.117003514, 0.117003514},
 {0.193716983, 0.193716983, 0.161601202, 0.161601202, 0.193716983, 0.193716983},
 {0.2556529, 0.2556529, 0.230551593, 0.230551593, 0.2556529, 0.2556529},
 {0.307686066, 0.307686066, 0.28800865, 0.28800865, 0.307686066, 0.307686066},
 {0.35173019, 0.35173019, 0.336045787, 0.336045787, 0.35173019, 0.35173019},
 {0.389243417, 0.389243417, 0.376527564, 0.376527564, 0.389243417, 0.389243417},
 {0.421468337, 0.421468337, 0.410985259, 0.410985259, 0.421468337, 0.421468337},
 {0.449453624, 0.449453624, 0.440663694, 0.440663694, 0.449453624, 0.449453624},
 {0.474072401, 0.474072401, 0.46657284, 0.46657284, 0.474072401, 0.474072401},
 {0.496005712, 0.496005712, 0.489499281, 0.489499281, 0.496005712, 0.496005712},
 {0.515769612, 0.515769612, 0.510038748, 0.510038748, 0.515769612, 0.515769612},
 {0.53380494, 0.53380494, 0.528680643, 0.528680643, 0.53380494, 0.53380494},
 {0.550471697, 0.550471697, 0.545822126, 0.545822126, 0.550471697, 0.550471697},
 {0.566063373, 0.566063373, 0.561785197, 0.561785197, 0.566063373, 0.566063373},
 {0.580818921, 0.580818921, 0.576830843, 0.576830843, 0.580818921, 0.580818921},
 {0.594932792, 0.594932792, 0.591170641, 0.591170641, 0.594932792, 0.594932792},
 {0.608563277, 0.608563277, 0.604976239, 0.604976239, 0.608563277, 0.608563277},
 {0.621839378, 0.621839378, 0.618387096, 0.618387096, 0.621839378, 0.621839378},
 {0.634955263, 0.634955263, 0.631591371, 0.631591371, 0.634955263, 0.634955263},
 {0.64912601, 0.64912601, 0.645584776, 0.645584776, 0.64912601, 0.64912601},
 {0.66752596, 0.66752596, 0.663086599, 0.663086599, 0.66752596, 0.66752596},
 {0.688594349, 0.688594349, 0.683556527, 0.683556527, 0.688594349, 0.688594349},
 {0.711680352, 0.711680352, 0.706235373, 0.706235373, 0.711680352, 0.711680352},
 {0.742115719, 0.742115719, 0.735177901, 0.735177901, 0.742115719, 0.742115719}};

    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);
}
