#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is simple two elements multi-domain example without boundary conditions. Initial
/// temperature and moisture distribution is not same in every node. This case should prove
/// that domain will try to reach equilibrium
////////////////////////////////////////////////////////////////////////////////////////////////////

class MultiDomain_2D_1 : public testing::Test
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

TEST_F(MultiDomain_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 20;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;
    const auto liquidPercent = 1.0;

    auto state = State(initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
    size_t nodeIndex = 0;
    auto T = 0.0;
    auto deltaT = 10.0;
    auto H = 0.0;
    auto deltaH = 0.1;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(
          nodeIndex,
          val,
          0.00,
          State(initialTemperature + T, initialMoistureContent + H, initialPressure, 0));
        ++nodeIndex;
        NodePool::Instance().createNode(
          nodeIndex,
          val,
          0.05,
          State(initialTemperature + T, initialMoistureContent + H, initialPressure, 0));
        T += deltaT;
        H += deltaH;
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

    const auto dTime = 360;
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
      {0.000944, 0.000944, 1.060053, 1.060053, 2.117104, 2.117104},
      {0.001987, 0.001987, 1.060083, 1.060083, 2.114588, 2.114588},
      {0.003113, 0.003113, 1.060099, 1.060099, 2.112348, 2.112348},
      {0.004305, 0.004305, 1.060107, 1.060107, 2.110308, 2.110308},
      {0.005551, 0.005551, 1.060109, 1.060109, 2.108417, 2.108417},
      {0.006839, 0.006839, 1.060109, 1.060109, 2.106638, 2.106638},
      {0.008159, 0.008159, 1.060107, 1.060107, 2.104943, 2.104943},
      {0.009505, 0.009505, 1.060103, 1.060103, 2.103312, 2.103312},
      {0.01087, 0.01087, 1.060099, 1.060099, 2.101731, 2.101731},
      {0.01225, 0.01225, 1.060094, 1.060094, 2.100187, 2.100187}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {22.291788, 22.291788, 30.001056, 30.001056, 37.718917, 37.718917},
      {24.058692, 24.058692, 30.002675, 30.002675, 35.95865, 35.95865},
      {25.42104, 25.42104, 30.004424, 30.004424, 34.600407, 34.600407},
      {26.471533, 26.471533, 30.006087, 30.006087, 33.552446, 33.552446},
      {27.281603, 27.281603, 30.007564, 30.007564, 32.743932, 32.743932},
      {27.906303, 27.906303, 30.008826, 30.008826, 32.120183, 32.120183},
      {28.388071, 28.388071, 30.009874, 30.009874, 31.638994, 31.638994},
      {28.75962, 28.75962, 30.01073, 30.01073, 31.267794, 31.267794},
      {29.046173, 29.046173, 30.01142, 30.01142, 30.98145, 30.98145},
      {29.267179, 29.267179, 30.01197, 30.01197, 30.760566, 30.760566}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
