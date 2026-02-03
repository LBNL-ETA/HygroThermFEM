#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ExcludeLiquidTransport_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        const auto relaxationParameter{0.8};
        const auto errorTolerance{1e-5};
        const auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ExcludeLiquidTransport_1, TestExample_1)
{
    const auto excludeWaterLiquidTransportation{true};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{false};
    const auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.99;
    const double initialPressure = 101325;
    constexpr auto liquidPercent = 1.0;

    auto state = HygroThermFEM::State(
      initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    // Material Properties (Cottaer Sandstone)
    constexpr double thermalConductivityDry{1.8};
    constexpr double density{2050.0};
    constexpr double porosity{0.22};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
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
      MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone",
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

    HygroThermFEM::MultiDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    constexpr auto hc = 5.0;
    constexpr auto airTemperature = 10.0;
    constexpr auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    domain.createBC_FixedHc(1, 2, bcCoeff);
    domain.createBC_FixedHc(5, 6, bcCoeff);

    constexpr auto dTime = 3600;
    constexpr auto nSteps = 24;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {62.287279, 62.287279, 63.000008, 63.000008, 62.287279, 62.287279},
      {61.553718, 61.553718, 62.999999, 62.999999, 61.553718, 61.553718},
      {60.798482, 60.798482, 62.99999, 62.99999, 60.798482, 60.798482},
      {60.020316, 60.020316, 62.999976, 62.999976, 60.020316, 60.020316},
      {59.217628, 59.217628, 62.999955, 62.999955, 59.217628, 59.217628},
      {58.38846, 58.38846, 62.999927, 62.999927, 58.38846, 58.38846},
      {57.530429, 57.530429, 62.999893, 62.999893, 57.530429, 57.530429},
      {56.640632, 56.640632, 62.999853, 62.999853, 56.640632, 56.640632},
      {55.715526, 55.715526, 62.999805, 62.999805, 55.715526, 55.715526},
      {54.750735, 54.750735, 62.99975, 62.99975, 54.750735, 54.750735},
      {53.740791, 53.740791, 62.999687, 62.999687, 53.740791, 53.740791},
      {52.678735, 52.678735, 62.999616, 62.999616, 52.678735, 52.678735},
      {51.555499, 51.555499, 62.999536, 62.999536, 51.555499, 51.555499},
      {50.35891, 50.35891, 62.999448, 62.999448, 50.35891, 50.35891},
      {49.071963, 49.071963, 62.999348, 62.999348, 49.071963, 49.071963},
      {47.669608, 47.669608, 62.999238, 62.999238, 47.669608, 47.669608},
      {46.112031, 46.112031, 62.999114, 62.999114, 46.112031, 46.112031},
      {44.328172, 44.328172, 62.998976, 62.998976, 44.328172, 44.328172},
      {42.461721, 42.461721, 62.998823, 62.998823, 42.461721, 42.461721},
      {40.555283, 40.555283, 62.998655, 62.998655, 40.555283, 40.555283},
      {38.606036, 38.606036, 62.998472, 62.998472, 38.606036, 38.606036},
      {36.610803, 36.610803, 62.998273, 62.998273, 36.610803, 36.610803},
      {34.565999, 34.565999, 62.998059, 62.998059, 34.565999, 34.565999},
      {32.467541, 32.467541, 62.997827, 62.997827, 32.467541, 32.467541}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.104113, 0.104113, 0.071611, 0.071611, 0.104113, 0.104113},
      {0.175255, 0.175255, 0.142905, 0.142905, 0.175255, 0.175255},
      {0.231249, 0.231249, 0.203676, 0.203676, 0.231249, 0.231249},
      {0.277033, 0.277033, 0.254138, 0.254138, 0.277033, 0.277033},
      {0.314938, 0.314938, 0.29596, 0.29596, 0.314938, 0.314938},
      {0.346562, 0.346562, 0.330765, 0.330765, 0.346562, 0.346562},
      {0.373154, 0.373154, 0.359918, 0.359918, 0.373154, 0.373154},
      {0.395718, 0.395718, 0.384535, 0.384535, 0.395718, 0.395718},
      {0.415072, 0.415072, 0.405531, 0.405531, 0.415072, 0.415072},
      {0.431886, 0.431886, 0.423649, 0.423649, 0.431886, 0.431886},
      {0.446706, 0.446706, 0.439498, 0.439498, 0.446706, 0.446706},
      {0.459986, 0.459986, 0.453582, 0.453582, 0.459986, 0.459986},
      {0.472106, 0.472106, 0.466317, 0.466317, 0.472106, 0.472106},
      {0.483388, 0.483388, 0.478058, 0.478058, 0.483388, 0.483388},
      {0.494117, 0.494117, 0.48911, 0.48911, 0.494117, 0.494117},
      {0.50456, 0.50456, 0.499753, 0.499753, 0.50456, 0.50456},
      {0.514988, 0.514988, 0.510263, 0.510263, 0.514988, 0.514988},
      {0.525723, 0.525723, 0.520947, 0.520947, 0.525723, 0.525723},
      {0.536765, 0.536765, 0.531892, 0.531892, 0.536765, 0.536765},
      {0.548091, 0.548091, 0.543114, 0.543114, 0.548091, 0.548091},
      {0.559707, 0.559707, 0.554624, 0.554624, 0.559707, 0.559707},
      {0.571624, 0.571624, 0.566432, 0.566432, 0.571624, 0.571624},
      {0.583857, 0.583857, 0.578554, 0.578554, 0.583857, 0.583857},
      {0.596424, 0.596424, 0.591005, 0.591005, 0.596424, 0.596424}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
