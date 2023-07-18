#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_VariableTARPHc_MultiTimestepBC : public testing::Test
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

TEST_F(MultiDomain_2D_VariableTARPHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::TARPCoefficients> bcCoeff{{20.0, 0.6},
                                                               {20.0, 0.5},
                                                               {20.0, 0.4},
                                                               {20.0, 0.3},
                                                               {20.0, 0.2},
                                                               {18.0, 0.2},
                                                               {16.0, 0.2},
                                                               {14.0, 0.2},
                                                               {12.0, 0.2},
                                                               {10.0, 0.2}};

    const auto surfaceTilt{90.0};

    createBC_TARPHc(domain, 1, 2, bcCoeff, surfaceTilt);

    const auto dTime = 3600;
    const auto nSteps = 10;

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    HygroThermFEM::TransientSubstitutionSolver solver{domain};
    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {2.962591, 2.962591, 0.003664, 0.003664, 0.000009, 0.000009},
      {4.409299, 4.409299, 0.009581, 0.009581, 0.000035, 0.000035},
      {5.131025, 5.131025, 0.016977, 0.016977, 0.000085, 0.000085},
      {5.452173, 5.452173, 0.025202, 0.025202, 0.000165, 0.000165},
      {5.254478, 5.254478, 0.033638, 0.033638, 0.000278, 0.000278},
      {5.059500, 5.059500, 0.042063, 0.042063, 0.000425, 0.000425},
      {4.811524, 4.811524, 0.050289, 0.050289, 0.000605, 0.000605},
      {4.527274, 4.527274, 0.058162, 0.058162, 0.000817, 0.000817},
      {4.220785, 4.220785, 0.065565, 0.065565, 0.001058, 0.001058},
      {3.903781, 3.903781, 0.072418, 0.072418, 0.001325, 0.001325},
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
      {3.696622, 3.696622, 1.920131, 1.920131, 1.437081, 1.437081},
      {5.359277, 5.359277, 3.642653, 3.642653, 3.087792, 3.087792},
      {6.376815, 6.376815, 4.989258, 4.989258, 4.510897, 4.510897},
      {7.075910, 7.075910, 6.009417, 6.009417, 5.632421, 5.632421},
      {7.563317, 7.563317, 6.766152, 6.766152, 6.480923, 6.480923},
      {7.909440, 7.909440, 7.321657, 7.321657, 7.110136, 7.110136},
      {8.125347, 8.125347, 7.710471, 7.710471, 7.559425, 7.559425},
      {8.218634, 8.218634, 7.953762, 7.953762, 7.854537, 7.854537},
      {8.200384, 8.200384, 8.068054, 8.068054, 8.014318, 8.014318},
      {8.080523, 8.080523, 8.066742, 8.066742, 8.053532, 8.053532},
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
