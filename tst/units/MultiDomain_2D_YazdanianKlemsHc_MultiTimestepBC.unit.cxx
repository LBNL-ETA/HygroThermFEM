#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_YazdanianKlemsHc_MultiTimestepBC : public testing::Test
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

TEST_F(MultiDomain_2D_YazdanianKlemsHc_MultiTimestepBC, TestExample_1)
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

    using HygroThermFEM::WindDirection;

    // Variable boundary conditions (temperature, humidity, wind speed and wind direction) over ten
    // timesteps.
    const std::vector<HygroThermFEM::YazdanianKlemsCoefficients> bcCoeff{
      {20.0, 0.6, 3.0, WindDirection::Windward},
      {20.0, 0.5, 3.0, WindDirection::Windward},
      {20.0, 0.4, 3.0, WindDirection::Windward},
      {20.0, 0.3, 4.0, WindDirection::Windward},
      {20.0, 0.2, 4.2, WindDirection::Windward},
      {18.0, 0.2, 4.6, WindDirection::Leeward},
      {16.0, 0.2, 5.0, WindDirection::Leeward},
      {14.0, 0.2, 5.3, WindDirection::Leeward},
      {12.0, 0.2, 5.5, WindDirection::Leeward},
      {10.0, 0.2, 5.9, WindDirection::Leeward}};

    createBC_YazdanianKlemsHc(domain, 1, 2, bcCoeff);

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
      {4.119315, 4.119315, 0.005435, 0.005435, 0.000014, 0.000014},
      {6.250389, 6.250389, 0.014202, 0.014202, 0.000058, 0.000058},
      {7.302360, 7.302360, 0.024890, 0.024890, 0.000145, 0.000145},
      {6.784879, 6.784879, 0.036153, 0.036153, 0.000285, 0.000285},
      {5.107673, 5.107673, 0.046459, 0.046459, 0.000480, 0.000480},
      {4.352752, 4.352752, 0.055703, 0.055703, 0.000725, 0.000725},
      {3.673174, 3.673174, 0.063756, 0.063756, 0.001013, 0.001013},
      {3.107434, 3.107434, 0.070659, 0.070659, 0.001337, 0.001337},
      {2.663352, 2.663352, 0.076546, 0.076546, 0.001688, 0.001688},
      {2.333544, 2.333544, 0.081553, 0.081553, 0.002056, 0.002056},
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
      {6.290786, 6.290786, 3.267036, 3.267036, 2.445144, 2.445144},
      {8.607262, 8.607262, 5.931511, 5.931511, 5.054433, 5.054433},
      {9.762987, 9.762987, 7.804784, 7.804784, 7.112859, 7.112859},
      {10.832775, 10.832775, 9.285148, 9.285148, 8.738640, 8.738640},
      {11.490443, 11.490443, 10.357788, 10.357788, 9.950430, 9.950430},
      {11.749035, 11.749035, 11.026079, 11.026079, 10.755448, 10.755448},
      {11.894193, 11.894193, 11.440727, 11.440727, 11.268301, 11.268301},
      {11.837916, 11.837916, 11.623753, 11.623753, 11.534301, 11.534301},
      {11.567187, 11.567187, 11.582028, 11.582028, 11.569993, 11.569993},
      {11.075399, 11.075399, 11.316698, 11.316698, 11.380395, 11.380395},
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
