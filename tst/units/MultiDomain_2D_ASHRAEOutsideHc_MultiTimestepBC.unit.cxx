#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_ASHRAEOutsideHc_MultiTimestepBC : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(MultiDomain_2D_ASHRAEOutsideHc_MultiTimestepBC, TestExample_1)
{
    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.0;
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
    for(size_t i = 1; i <= (HygroThermFEM::maxNodeIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        createElement(domain, node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature, humidity and wind speed) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEOutsideCoefficients> bcCoeff{{20.0, 0.6, 3},
                                                                        {20.0, 0.5, 3},
                                                                        {20.0, 0.4, 3},
                                                                        {20.0, 0.3, 4},
                                                                        {20.0, 0.2, 4.2},
                                                                        {18.0, 0.2, 4.6},
                                                                        {16.0, 0.2, 5},
                                                                        {14.0, 0.2, 5.3},
                                                                        {12.0, 0.2, 5.5},
                                                                        {10.0, 0.2, 5.9}};

    createBC_ASHRAEOutsideHc(domain, 1, 2, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 10;

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        HygroThermFEM::TransientSubstitutionSolver solver;
        auto aSolution = solver.transient(domain, temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {5.839487, 5.839487, 0.008245, 0.008245, 0.000024, 0.000024},
      {8.360881, 8.360881, 0.020737, 0.020737, 0.000100, 0.000100},
      {7.723775, 7.723775, 0.034390, 0.034390, 0.000248, 0.000248},
      {5.210149, 5.210149, 0.046633, 0.046633, 0.000473, 0.000473},
      {3.487294, 3.487294, 0.055719, 0.055719, 0.000766, 0.000766},
      {2.673545, 2.673545, 0.063060, 0.063060, 0.001116, 0.001116},
      {2.299217, 2.299217, 0.069381, 0.069381, 0.001509, 0.001509},
      {2.098220, 2.098220, 0.074978, 0.074978, 0.001929, 0.001929},
      {1.975006, 1.975006, 0.079962, 0.079962, 0.002360, 0.002360},
      {1.895041, 1.895041, 0.084395, 0.084395, 0.002787, 0.002787},
    };

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {10.519157, 10.519157, 5.461415, 5.461415, 4.087477, 4.087477},
      {12.721219, 12.721219, 9.048708, 9.048708, 7.800587, 7.800587},
      {13.533831, 13.533831, 11.212281, 11.212281, 10.353968, 10.353968},
      {14.527596, 14.527596, 12.820346, 12.820346, 12.199840, 12.199840},
      {15.451696, 15.451696, 14.105054, 14.105054, 13.625716, 13.625716},
      {15.794810, 15.794810, 14.919423, 14.919423, 14.593922, 14.593922},
      {15.460194, 15.460194, 15.157125, 15.157125, 15.015402, 15.015402},
      {14.620576, 14.620576, 14.859222, 14.859222, 14.898481, 14.898481},
      {13.435956, 13.435956, 14.124345, 14.124345, 14.319073, 14.319073},
      {11.978456, 11.978456, 13.034403, 13.034403, 13.357576, 13.357576},
    };

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
