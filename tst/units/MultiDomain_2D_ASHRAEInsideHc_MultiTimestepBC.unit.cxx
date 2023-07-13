#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_ASHRAEInsideHc_MultiTimestepBC : public testing::Test
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

TEST_F(MultiDomain_2D_ASHRAEInsideHc_MultiTimestepBC, TestExample_1)
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

    // Variable boundary conditions (temperature, pressure and humidity) over ten timesteps.
    const std::vector<HygroThermFEM::ASHRAEInsideCoefficients> bcCoeff{{20.0, 0.6, 101325.0},
                                                                       {20.0, 0.5, 101325.0},
                                                                       {20.0, 0.4, 101325.0},
                                                                       {20.0, 0.3, 101325.0},
                                                                       {20.0, 0.2, 101325.0},
                                                                       {18.0, 0.2, 101325.0},
                                                                       {16.0, 0.2, 101325.0},
                                                                       {14.0, 0.2, 101325.0},
                                                                       {12.0, 0.2, 101325.0},
                                                                       {10.0, 0.2, 101325.0}};

    const auto surfaceTilt{90.0};    // degrees
    const auto surfaceHeight{1.0};   // meters

    createBC_ASHRAEInsideHc(domain, 1, 2, bcCoeff, surfaceHeight, surfaceTilt);

    const auto dTime = 3600;
    const auto nSteps = 10;

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    auto humidities = properties(HygroThermFEM::Variable::humidity);
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
      {2.779607, 2.779607, 0.003409, 0.003409, 0.000008, 0.000008},
      {4.180510, 4.180510, 0.008930, 0.008930, 0.000032, 0.000032},
      {4.907598, 4.907598, 0.015857, 0.015857, 0.000078, 0.000078},
      {5.187634, 5.187634, 0.023601, 0.023601, 0.000151, 0.000151},
      {5.116001, 5.116001, 0.031605, 0.031605, 0.000254, 0.000254},
      {4.977129, 4.977129, 0.039663, 0.039663, 0.000388, 0.000388},
      {4.794767, 4.794767, 0.047614, 0.047614, 0.000553, 0.000553},
      {4.588008, 4.588008, 0.055336, 0.055336, 0.000749, 0.000749},
      {4.373404, 4.373404, 0.062748, 0.062748, 0.000972, 0.000972},
      {4.168542, 4.168542, 0.069817, 0.069817, 0.001221, 0.001221},
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
      {3.341972, 3.341972, 1.735960, 1.735960, 1.299243, 1.299243},
      {4.910435, 4.910435, 3.327176, 3.327176, 2.817004, 2.817004},
      {5.900680, 5.900680, 4.596293, 4.596293, 4.148669, 4.148669},
      {6.596696, 6.596696, 5.575783, 5.575783, 5.216752, 5.216752},
      {7.091261, 7.091261, 6.314985, 6.314985, 6.038688, 6.038688},
      {7.419937, 7.419937, 6.851802, 6.851802, 6.647231, 6.647231},
      {7.611634, 7.611634, 7.218782, 7.218782, 7.074979, 7.074979},
      {7.689663, 7.689663, 7.443672, 7.443672, 7.350901, 7.350901},
      {7.675291, 7.675291, 7.551010, 7.551010, 7.500649, 7.500649},
      {7.591926, 7.591926, 7.564887, 7.564887, 7.548706, 7.548706},
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
