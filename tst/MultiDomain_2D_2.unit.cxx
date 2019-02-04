#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

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

    auto state =
      MoisThermFEM::State(initialTemperature, initialMoistureContent, initialPressure, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
    }

    auto & material = MaterialPool::Instance().createMaterial(
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

    MoisThermFEM::MultiDomain domain;

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
    const auto humidity = 0.2;

    domain.createMoistureBC(1, 2, hc, airTemperature, humidity);

    const auto dTime = 3600;
    const auto nSteps = 10;

    auto temperatures = NodePool::Instance().properties(MoisThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(MoisThermFEM::Variable::humidity);
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
      {1.3590234, 1.3590234, 0.0015533155, 0.0015533155, 3.4745566e-006, 3.4745566e-006},
      {2.455398, 2.455398, 0.0044425073, 0.0044425073, 1.366852e-005, 1.366852e-005},
      {3.3308922, 3.3308922, 0.0084617234, 0.0084617234, 3.3582494e-005, 3.3582494e-005},
      {4.0211014, 4.0211014, 0.013430135, 0.013430135, 6.5979379e-005, 6.5979379e-005},
      {4.5566542, 4.5566542, 0.019190124, 0.019190124, 0.00011339779, 0.00011339779},
      {4.9639115, 4.9639115, 0.025604801, 0.025604801, 0.00017816616, 0.00017816616},
      {5.265464, 5.265464, 0.032555882, 0.032555882, 0.00026241724, 0.00026241724},
      {5.4893015, 5.4893015, 0.039295983, 0.039295983, 0.00036638299, 0.00036638299},
      {5.6825361, 5.6825361, 0.046298788, 0.046298788, 0.00049154203, 0.00049154203},
      {5.8199264, 5.8199264, 0.053530535, 0.053530535, 0.00063931787, 0.00063931787}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.65796851, 0.65796851, 0.34200269, 0.34200269, 0.25596504, 0.25596504},
      {1.1011411, 1.1011411, 0.72533834, 0.72533834, 0.60725804, 0.60725804},
      {1.4906398, 1.4906398, 1.1076894, 1.1076894, 0.98179576, 0.98179576},
      {1.8619457, 1.8619457, 1.4832841, 1.4832841, 1.3571245, 1.3571245},
      {2.2233285, 2.2233285, 1.8514668, 1.8514668, 1.7271049, 1.7271049},
      {2.5768635, 2.5768635, 2.2122834, 2.2122834, 2.0902268, 2.0902268},
      {2.9231654, 2.9231654, 2.565871, 2.565871, 2.4462129, 2.4462129},
      {3.2625072, 3.2625072, 2.9123797, 2.9123797, 2.7951058, 2.7951058},
      {3.5950499, 3.5950499, 3.2519536, 3.2519536, 3.137024, 3.137024},
      {3.9209386, 3.9209386, 3.5847337, 3.5847337, 3.4721029, 3.4721029}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
