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
    const auto liquidPercent = 1.0;

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
    const auto hc = 5.0;
    const auto airTemperature = 10.0;
    const auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    domain.createMoistureBCFixedHc(1, 2, bcCoeff);
    domain.createMoistureBCFixedHc(5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 24;

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
      {62.287230, 62.287230, 63.000009, 63.000009, 62.287230, 62.287230},
      {61.553637, 61.553637, 62.999999, 62.999999, 61.553637, 61.553637},
      {60.798375, 60.798375, 62.999991, 62.999991, 60.798375, 60.798375},
      {60.020188, 60.020188, 62.999976, 62.999976, 60.020188, 60.020188},
      {59.217482, 59.217482, 62.999955, 62.999955, 59.217482, 59.217482},
      {58.388298, 58.388298, 62.999927, 62.999927, 58.388298, 58.388298},
      {57.530252, 57.530252, 62.999893, 62.999893, 57.530252, 57.530252},
      {56.640443, 56.640443, 62.999853, 62.999853, 56.640443, 56.640443},
      {55.715322, 55.715322, 62.999805, 62.999805, 55.715322, 55.715322},
      {54.750518, 54.750518, 62.999750, 62.999750, 54.750518, 54.750518},
      {53.740559, 53.740559, 62.999687, 62.999687, 53.740559, 53.740559},
      {52.678487, 52.678487, 62.999616, 62.999616, 52.678487, 52.678487},
      {51.555232, 51.555232, 62.999537, 62.999537, 51.555232, 51.555232},
      {50.358620, 50.358620, 62.999448, 62.999448, 50.358620, 50.358620},
      {49.071646, 49.071646, 62.999349, 62.999349, 49.071646, 49.071646},
      {47.669254, 47.669254, 62.999238, 62.999238, 47.669254, 47.669254},
      {46.111628, 46.111628, 62.999115, 62.999115, 46.111628, 46.111628},
      {44.327694, 44.327694, 62.998976, 62.998976, 44.327694, 44.327694},
      {42.461224, 42.461224, 62.998823, 62.998823, 42.461224, 42.461224},
      {40.554766, 40.554766, 62.998655, 62.998655, 40.554766, 40.554766},
      {38.605496, 38.605496, 62.998472, 62.998472, 38.605496, 38.605496},
      {36.610240, 36.610240, 62.998274, 62.998274, 36.610240, 36.610240},
      {34.565409, 34.565409, 62.998059, 62.998059, 34.565409, 34.565409},
      {32.466922, 32.466922, 62.997828, 62.997828, 32.466922, 32.466922}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.104018, 0.104018, 0.071546, 0.071546, 0.104018, 0.104018},
      {0.175126, 0.175126, 0.142796, 0.142796, 0.175126, 0.175126},
      {0.231103, 0.231103, 0.203542, 0.203542, 0.231103, 0.231103},
      {0.276881, 0.276881, 0.253991, 0.253991, 0.276881, 0.276881},
      {0.314787, 0.314787, 0.295811, 0.295811, 0.314787, 0.314787},
      {0.346417, 0.346417, 0.330618, 0.330618, 0.346417, 0.346417},
      {0.373016, 0.373016, 0.359777, 0.359777, 0.373016, 0.373016},
      {0.395589, 0.395589, 0.384403, 0.384403, 0.395589, 0.395589},
      {0.414954, 0.414954, 0.405408, 0.405408, 0.414954, 0.414954},
      {0.431777, 0.431777, 0.423535, 0.423535, 0.431777, 0.431777},
      {0.446607, 0.446607, 0.439395, 0.439395, 0.446607, 0.446607},
      {0.459896, 0.459896, 0.453487, 0.453487, 0.459896, 0.459896},
      {0.472024, 0.472024, 0.466231, 0.466231, 0.472024, 0.472024},
      {0.483312, 0.483312, 0.477979, 0.477979, 0.483312, 0.483312},
      {0.494048, 0.494048, 0.489038, 0.489038, 0.494048, 0.494048},
      {0.504496, 0.504496, 0.499687, 0.499687, 0.504496, 0.504496},
      {0.514928, 0.514928, 0.510201, 0.510201, 0.514928, 0.514928},
      {0.525666, 0.525666, 0.520889, 0.520889, 0.525666, 0.525666},
      {0.536710, 0.536710, 0.531837, 0.531837, 0.536710, 0.536710},
      {0.548038, 0.548038, 0.543061, 0.543061, 0.548038, 0.548038},
      {0.559656, 0.559656, 0.554572, 0.554572, 0.559656, 0.559656},
      {0.571573, 0.571573, 0.566381, 0.566381, 0.571573, 0.571573},
      {0.583806, 0.583806, 0.578502, 0.578502, 0.583806, 0.583806},
      {0.596373, 0.596373, 0.590953, 0.590953, 0.596373, 0.596373}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
