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

    SimulationProperties::Instance().setCalculationParameters(excludeWaterLiquidTransportation,
                                                              excludeHeatOfEvaporation,
                                                              excludeCapillaryConduction,
                                                              excludeVaporDiffusionConduction);

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

    auto & material = MaterialPool::Instance().createSolidMaterial(
      "Cottaer Sandstone",
      2050,                       /// density
      0.22,                       /// porosity
      850,                        /// specific heat capacity (dry)
      15,                         /// diffusion resistance factor
      {{0.0, 1.8}, {180, 1.8}},   /// thermal conductivity as function of water content
      {{0, 0},                    /// liquid transportation coefficient
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// sorption curve
       {0.5, 5.3},
       {0.65, 8.4},
       {0.8, 12},
       {0.93, 17},
       {0.95, 25},
       {0.99, 63},
       {0.995, 83},
       {0.999, 120},
       {1, 180}});

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

    domain.createMoistureBCFixedHc(1, 2, airTemperature, hc, humidity);
    domain.createMoistureBCFixedHc(5, 6, airTemperature, hc, humidity);

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
      {61.553636, 61.553636, 62.999999, 62.999999, 61.553636, 61.553636},
      {60.798373, 60.798373, 62.999991, 62.999991, 60.798373, 60.798373},
      {60.020185, 60.020185, 62.999976, 62.999976, 60.020185, 60.020185},
      {59.217480, 59.217480, 62.999955, 62.999955, 59.217480, 59.217480},
      {58.388297, 58.388297, 62.999927, 62.999927, 58.388297, 58.388297},
      {57.530255, 57.530255, 62.999893, 62.999893, 57.530255, 57.530255},
      {56.640453, 56.640453, 62.999853, 62.999853, 56.640453, 56.640453},
      {55.715346, 55.715346, 62.999805, 62.999805, 55.715346, 55.715346},
      {54.750562, 54.750562, 62.999750, 62.999750, 54.750562, 54.750562},
      {53.740636, 53.740636, 62.999687, 62.999687, 53.740636, 53.740636},
      {52.678612, 52.678612, 62.999616, 62.999616, 52.678612, 52.678612},
      {51.555426, 51.555426, 62.999537, 62.999537, 51.555426, 51.555426},
      {50.358912, 50.358912, 62.999448, 62.999448, 50.358912, 50.358912},
      {49.072078, 49.072078, 62.999349, 62.999349, 49.072078, 49.072078},
      {47.669890, 47.669890, 62.999238, 62.999238, 47.669890, 47.669890},
      {46.112571, 46.112571, 62.999115, 62.999115, 46.112571, 46.112571},
      {44.329132, 44.329132, 62.998977, 62.998977, 44.329132, 44.329132},
      {42.463062, 42.463062, 62.998824, 62.998824, 42.463062, 42.463062},
      {40.557090, 40.557090, 62.998656, 62.998656, 40.557090, 40.557090},
      {38.608405, 38.608405, 62.998473, 62.998473, 38.608405, 38.608405},
      {36.613843, 36.613843, 62.998275, 62.998275, 36.613843, 36.613843},
      {34.569835, 34.569835, 62.998060, 62.998060, 34.569835, 34.569835},
      {32.472316, 32.472316, 62.997829, 62.997829, 32.472316, 32.472316}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.104023, 0.104023, 0.071551, 0.071551, 0.104023, 0.104023},
      {0.175137, 0.175137, 0.142808, 0.142808, 0.175137, 0.175137},
      {0.231119, 0.231119, 0.203559, 0.203559, 0.231119, 0.231119},
      {0.276894, 0.276894, 0.254006, 0.254006, 0.276894, 0.276894},
      {0.314788, 0.314788, 0.295813, 0.295813, 0.314788, 0.314788},
      {0.346393, 0.346393, 0.330593, 0.330593, 0.346393, 0.346393},
      {0.372953, 0.372953, 0.359708, 0.359708, 0.372953, 0.372953},
      {0.395469, 0.395469, 0.384272, 0.384272, 0.395469, 0.395469},
      {0.414758, 0.414758, 0.405193, 0.405193, 0.414758, 0.414758},
      {0.431483, 0.431483, 0.423211, 0.423211, 0.431483, 0.431483},
      {0.446190, 0.446190, 0.438935, 0.438935, 0.446190, 0.446190},
      {0.459329, 0.459329, 0.452861, 0.452861, 0.459329, 0.459329},
      {0.471275, 0.471275, 0.465404, 0.465404, 0.471275, 0.471275},
      {0.482347, 0.482347, 0.476912, 0.476912, 0.482347, 0.482347},
      {0.492824, 0.492824, 0.487685, 0.487685, 0.492824, 0.492824},
      {0.502964, 0.502964, 0.497993, 0.497993, 0.502964, 0.502964},
      {0.513026, 0.513026, 0.508098, 0.508098, 0.513026, 0.513026},
      {0.523315, 0.523315, 0.518290, 0.518290, 0.523315, 0.523315},
      {0.533874, 0.533874, 0.528701, 0.528701, 0.533874, 0.533874},
      {0.544689, 0.544689, 0.539355, 0.539355, 0.544689, 0.544689},
      {0.555760, 0.555760, 0.550261, 0.550261, 0.555760, 0.555760},
      {0.567094, 0.567094, 0.561424, 0.561424, 0.567094, 0.567094},
      {0.578703, 0.578703, 0.572854, 0.572854, 0.578703, 0.578703},
      {0.590600, 0.590600, 0.584565, 0.584565, 0.590600, 0.590600}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
