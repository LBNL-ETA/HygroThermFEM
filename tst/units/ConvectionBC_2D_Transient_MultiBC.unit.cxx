#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class ConvectionBC_2D_Transient_MultiBC : public testing::Test
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

TEST_F(ConvectionBC_2D_Transient_MultiBC, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Three elements with simple convection BC.");

    NodePool::Instance().createNode(1, 0.2, 0.05);
    NodePool::Instance().createNode(2, 0.2, 0.00);
    NodePool::Instance().createNode(3, 0.1, 0.05);
    NodePool::Instance().createNode(4, 0.1, 0.00);
    NodePool::Instance().createNode(5, 0.0, 0.05);
    NodePool::Instance().createNode(6, 0.0, 0.00);

    // Material Properties
    const double thermalConductivityDry{1.0};
    const double density{2050.0};
    const double porosity{0.0};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.0}, {180, 1.0}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.0}, {1, 1.0}};
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
      MaterialPool::Instance().createSolidMaterial("Test Material",
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

    HygroThermFEM::ThermalDomain domain;

    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions

    // Here we create boundary conditions for every time-step
    const std::vector<HygroThermFEM::FixedBCHCCoefficients> bcCoeffsTransient{
      {20.0, 2.4}, {20.0, 2.2}, {20.0, 2.0}, {18.0, 1.5}};
    //{20.0, 2.4}, {20.0, 2.2}, {20.0, 2.0}};

    const auto hc2 = 15.0;
    const auto temperatureAir2 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2};

    domain.createBC_FixedHc(1, 2, bcCoeffsTransient);
    domain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto dTime = 3600;
    const auto nSteps = 4;

    auto temperatures = properties(HygroThermFEM::Variable::temperature);
    std::vector<std::vector<double>> temperaturesSolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;
    size_t timestepIndex{0};

    for(unsigned i = 0; i < nSteps; ++i)
    {
        temperatures = transient(domain, temperatures, dTime, timestepIndex).solution;
        temperaturesSolution.push_back(temperatures);
        fluxSolution.push_back(domain.flux());
        ++timestepIndex;
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.132130505, 1.132130505, -0.656245111, -0.656245111, -5.621029352, -5.621029352},
      {1.552647394, 1.552647394, -1.488060905, -1.488060905, -8.554988708, -8.554988708},
      {1.502118281, 1.502118281, -2.319745532, -2.319745532, -10.16719396, -10.16719396},
      {0.895825525, 0.895825525, -3.137113323, -3.137113323, -11.12633933, -11.12633933}};

    std::vector<std::vector<HygroThermFEM::NodeFlux>> correctFluxSolution{{{-17.88375616, 0},
                                                                           {-17.88375616, 0},
                                                                           {-33.76579929, 0},
                                                                           {-33.76579929, 0},
                                                                           {-49.64784241, 0},
                                                                           {-49.64784241, 0}},
                                                                          {{-30.40708298, 0},
                                                                           {-30.40708298, 0},
                                                                           {-50.53818051, 0},
                                                                           {-50.53818051, 0},
                                                                           {-70.66927803, 0},
                                                                           {-70.66927803, 0}},
                                                                          {{-38.21863813, 0},
                                                                           {-38.21863813, 0},
                                                                           {-58.34656121, 0},
                                                                           {-58.34656121, 0},
                                                                           {-78.47448429, 0},
                                                                           {-78.47448429, 0}},
                                                                          {{-40.32938848, 0},
                                                                           {-40.32938848, 0},
                                                                           {-60.11082427, 0},
                                                                           {-60.11082427, 0},
                                                                           {-79.89226005, 0},
                                                                           {-79.89226005, 0}}};

    EXPECT_EQ(temperaturesSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperaturesSolution[i][j], 1e-6);
            EXPECT_NEAR(correctFluxSolution[i][j].x, fluxSolution[i][j].x, 1e-6);
            EXPECT_NEAR(correctFluxSolution[i][j].y, fluxSolution[i][j].y, 1e-6);
        }
    }
}
