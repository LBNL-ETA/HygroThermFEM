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
    {}

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
    const double initialMoistureContent = 0.8;
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
    const auto airTemperature = 20.0;
    const auto humidity = 0.0;

    domain.createMoistureBCFixedHc(1, 2, airTemperature, hc, humidity);

    const auto dTime = 3600;
    const auto nSteps = 10;

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
      {9.926993, 9.926993, 11.999349, 11.999349, 12.000188, 12.000188},
      {7.912585, 7.912585, 11.997777, 11.997777, 12.000457, 12.000457},
      {6.117678, 6.117678, 11.995133, 11.995133, 12.000752, 12.000752},
      {4.805022, 4.805022, 11.991280, 11.991280, 12.001058, 12.001058},
      {3.917136, 3.917136, 11.986095, 11.986095, 12.001369, 12.001369},
      {3.080759, 3.080759, 11.979486, 11.979486, 12.001682, 12.001682},
      {2.357192, 2.357192, 11.971466, 11.971466, 12.001994, 12.001994},
      {1.783685, 1.783685, 11.962139, 11.962139, 12.002297, 12.002297},
      {1.336457, 1.336457, 11.951589, 11.951589, 12.002585, 12.002585},
      {0.993381, 0.993381, 11.939897, 11.939897, 12.002847, 12.002847}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.753716, 1.753716, 0.886138, 0.886138, 0.652476, 0.652476},
      {2.857198, 2.857198, 1.851376, 1.851376, 1.535243, 1.535243},
      {3.810069, 3.810069, 2.799775, 2.799775, 2.466334, 2.466334},
      {4.724913, 4.724913, 3.729253, 3.729253, 3.396236, 3.396236},
      {5.630198, 5.630198, 4.646860, 4.646860, 4.317084, 4.317084},
      {6.533387, 6.533387, 5.557918, 5.557918, 5.230723, 5.230723},
      {7.426514, 7.426514, 6.460546, 6.460546, 6.136255, 6.136255},
      {8.295243, 8.295243, 7.346675, 7.346675, 7.027501, 7.027501},
      {9.131525, 9.131525, 8.208501, 8.208501, 7.897086, 7.897086},
      {9.929798, 9.929798, 9.039412, 9.039412, 8.738197, 8.738197}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
