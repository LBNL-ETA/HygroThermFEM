#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_VariableHc_MultiTimestepBC : public testing::Test
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

TEST_F(MultiDomain_2D_VariableHc_MultiTimestepBC, TestExample_1)
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
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions

    // Variable boundary conditions (temperature and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::VariableBCHCCoefficients> bcCoeff{{20.0, 0.6},
                                                                       {20.0, 0.5},
                                                                       {20.0, 0.4},
                                                                       {20.0, 0.3},
                                                                       {20.0, 0.2},
                                                                       {18.0, 0.2},
                                                                       {16.0, 0.2},
                                                                       {14.0, 0.2},
                                                                       {12.0, 0.2},
                                                                       {10.0, 0.2}};

    domain.createMoistureBCVariableHc(1, 2, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 10;

    auto temperatures = NodePool::Instance().properties(HygroThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;
    size_t timestepIndex{0u};

    for(auto i = 0; i < nSteps; ++i)
    {
        auto aSolution = domain.transient(temperatures, humidities, dTime, timestepIndex);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {2.962928, 2.962928, 0.003664, 0.003664, 8e-06, 8e-06},
      {4.409729, 4.409729, 0.009582, 0.009582, 3.4e-05, 3.4e-05},
      {5.131445, 5.131445, 0.016979, 0.016979, 8.4e-05, 8.4e-05},
      {5.452859, 5.452859, 0.025204, 0.025204, 0.000164, 0.000164},
      {5.254719, 5.254719, 0.033641, 0.033641, 0.000277, 0.000277},
      {5.059675, 5.059675, 0.042067, 0.042067, 0.000424, 0.000424},
      {4.811645, 4.811645, 0.050293, 0.050293, 0.000605, 0.000605},
      {4.527345, 4.527345, 0.058167, 0.058167, 0.000817, 0.000817},
      {4.220815, 4.220815, 0.065571, 0.065571, 0.001058, 0.001058},
      {3.903787, 3.903787, 0.072424, 0.072424, 0.001324, 0.001324}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {3.697291, 3.697291, 1.920478, 1.920478, 1.437341, 1.437341},
      {5.360145, 5.360145, 3.643259, 3.643259, 3.088310, 3.088310},
      {6.377747, 6.377747, 4.990022, 4.990022, 4.511598, 4.511598},
      {7.076860, 7.076860, 6.010269, 6.010269, 5.633235, 5.633235},
      {7.564263, 7.564263, 6.767048, 6.767048, 6.481799, 6.481799},
      {7.910212, 7.910212, 7.322486, 7.322486, 7.110977, 7.110977},
      {8.126035, 8.126035, 7.711228, 7.711228, 7.560203, 7.560203},
      {8.219259, 8.219259, 7.954453, 7.954453, 7.855250, 7.855250},
      {8.200960, 8.200960, 8.068688, 8.068688, 8.014971, 8.014971},
      {8.081063, 8.081063, 8.067330, 8.067330, 8.054137, 8.054137}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
