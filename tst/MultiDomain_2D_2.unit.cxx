#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_2 : public testing::Test
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

TEST_F(MultiDomain_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

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
    const auto hc = 1.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.6;

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
      {1.0892812, 1.0892812, 0.0012989977, 0.0012989977, 2.9706892e-06, 2.9706892e-06},
      {2.112231, 2.112231, 0.0039597038, 0.0039597038, 1.2478773e-05, 1.2478773e-05},
      {3.0004066, 3.0004066, 0.0078840697, 0.0078840697, 3.2317775e-05, 3.2317775e-05},
      {3.7243508, 3.7243508, 0.01288011, 0.01288011, 6.6196586e-05, 6.6196586e-05},
      {4.3366186, 4.3366186, 0.018846943, 0.018846943, 0.00011788937, 0.00011788937},
      {4.8661864, 4.8661864, 0.025712567, 0.025712567, 0.00019125949, 0.00019125949},
      {5.3604117, 5.3604117, 0.033420827, 0.033420827, 0.00029025975, 0.00029025975},
      {6.1637797, 6.1637797, 0.041924682, 0.041924682, 0.00041891844, 0.00041891844},
      {6.8813336, 6.8813336, 0.051184035, 0.051184035, 0.00058132574, 0.00058132574},
      {7.5319337, 7.5319337, 0.061185046, 0.061185046, 0.00078168635, 0.00078168635}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.2789614, 1.2789614, 0.66454142, 0.66454142, 0.49736231, 0.49736231},
      {2.1030717, 2.1030717, 1.3898868, 1.3898868, 1.165353, 1.165353},
      {2.7994741, 2.7994741, 2.0924536, 2.0924536, 1.85922, 1.85922},
      {3.4431536, 3.4431536, 2.7630832, 2.7630832, 2.5356936, 2.5356936},
      {4.0522411, 4.0522411, 3.4023248, 3.4023248, 3.1842994, 3.1842994},
      {4.6319356, 4.6319356, 4.0116952, 4.0116952, 3.8035377, 3.8035377},
      {5.1843204, 5.1843204, 4.5926044, 4.5926044, 4.3940863, 4.3940863},
      {5.7103707, 5.7103707, 5.1461109, 5.1461109, 4.9569083, 4.9569083},
      {6.2113848, 6.2113848, 5.6734198, 5.6734198, 5.4931476, 5.4931476},
      {6.6883167, 6.6883167, 6.1755938, 6.1755938, 6.0038875, 6.0038875}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
