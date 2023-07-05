#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class Moisture_2D_TwoElements_3 : public testing::Test
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

TEST_F(Moisture_2D_TwoElements_3, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.0;
    const auto initialPressure = 101325.0;
    const auto liquidPercent = 1.0;

    const State state(initialTemperature, initialHumidity, initialPressure, liquidPercent);

    NodePool::Instance().createNode(1, 0.15, 0.05, state);
    NodePool::Instance().createNode(2, 0.15, 0, state);
    NodePool::Instance().createNode(3, 0.05, 0.05, state);
    NodePool::Instance().createNode(4, 0.05, 0, state);
    NodePool::Instance().createNode(5, 0, 0.05, state);
    NodePool::Instance().createNode(6, 0, 0, state);

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

    HygroThermFEM::MoistureDomain domain;

    /// Create elements
    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto airTemperature = 20.0;
    const auto airHumidity = 0.5;

    const HygroThermFEM::TemperatureAndHumidity bcCoeff{airTemperature, airHumidity};

    domain.createBC_FixedHumidity(5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 24;

    auto humidities = properties(HygroThermFEM::Variable::humidity);
    std::vector<double> timesteps;
    std::vector<std::vector<double>> humiditySolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        auto solution = transient(domain, humidities, dTime);
        humidities = solution.solution;
        timesteps.push_back(solution.dTime);
        auto humidityContent = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
        humiditySolution.push_back(humidityContent);
        fluxSolution.push_back(domain.flux());
    }

    std::vector<double> correctTimesteps{3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600,
                                         3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600};

    std::vector<std::vector<double>> correctHumiditySolution{
      {2e-06, 2e-06, 0.001245, 0.001245, 0.5, 0.5},
      {7e-06, 7e-06, 0.002486, 0.002486, 0.5, 0.5},
      {1.4e-05, 1.4e-05, 0.003722, 0.003722, 0.5, 0.5},
      {2.4e-05, 2.4e-05, 0.004953, 0.004953, 0.5, 0.5},
      {3.6e-05, 3.6e-05, 0.00618, 0.00618, 0.5, 0.5},
      {5e-05, 5e-05, 0.007402, 0.007402, 0.5, 0.5},
      {6.7e-05, 6.7e-05, 0.00862, 0.00862, 0.5, 0.5},
      {8.6e-05, 8.6e-05, 0.009833, 0.009833, 0.5, 0.5},
      {0.000107, 0.000107, 0.011042, 0.011042, 0.5, 0.5},
      {0.000131, 0.000131, 0.012246, 0.012246, 0.5, 0.5},
      {0.000157, 0.000157, 0.013446, 0.013446, 0.5, 0.5},
      {0.000185, 0.000185, 0.014642, 0.014642, 0.5, 0.5},
      {0.000215, 0.000215, 0.015833, 0.015833, 0.5, 0.5},
      {0.000248, 0.000248, 0.017019, 0.017019, 0.5, 0.5},
      {0.000283, 0.000283, 0.018201, 0.018201, 0.5, 0.5},
      {0.00032, 0.00032, 0.019379, 0.019379, 0.5, 0.5},
      {0.00036, 0.00036, 0.020552, 0.020552, 0.5, 0.5},
      {0.000401, 0.000401, 0.021721, 0.021721, 0.5, 0.5},
      {0.000445, 0.000445, 0.022886, 0.022886, 0.5, 0.5},
      {0.000491, 0.000491, 0.024047, 0.024047, 0.5, 0.5},
      {0.000539, 0.000539, 0.025203, 0.025203, 0.5, 0.5},
      {0.000589, 0.000589, 0.026355, 0.026355, 0.5, 0.5},
      {0.000642, 0.000642, 0.027503, 0.027503, 0.5, 0.5},
      {0.000696, 0.000696, 0.028646, 0.028646, 0.5, 0.5}};


    EXPECT_EQ(correctTimesteps.size(), timesteps.size());

    for(size_t i = 0u; i < correctTimesteps.size(); ++i)
    {
        EXPECT_NEAR(correctTimesteps[i], timesteps[i], 1e-6);
    }

    EXPECT_EQ(correctHumiditySolution.size(), humiditySolution.size());

    for(auto i = 0u; i < humiditySolution.size(); ++i)
    {
        for(auto j = 0u; j < humiditySolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctHumiditySolution[i][j], humiditySolution[i][j], 1e-6);
        }
    }
}
