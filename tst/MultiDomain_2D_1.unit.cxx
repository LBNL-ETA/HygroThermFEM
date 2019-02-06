#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

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

    auto state = State(initialTemperature, initialMoistureContent, initialPressure, 0);
    size_t nodeIndex = 0;
    auto T = 0.0;
    auto deltaT = 10.0;
    //	auto deltaT = 10.0;
    auto H = 0.0;
    auto deltaH = 0.1;
    //	auto deltaH = 0.1;
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
    /// const auto dTime = 3600;
    const auto nSteps = 10;
    /// const auto nSteps = 10;

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
      {0.001452108, 0.001452108, 1.0607746, 1.0607746, 2.1163846, 2.1163846},
      {0.0029036488, 0.0029036488, 1.061336, 1.061336, 2.1133413, 2.1133413},
      {0.0043541072, 0.0043541072, 1.061747, 1.061747, 2.1107096, 2.1107096},
      {0.0058031169, 0.0058031169, 1.06205, 1.06205, 2.1083787, 2.1083787},
      {0.0072504252, 0.0072504252, 1.0622746, 1.0622746, 2.1062703, 2.1062703},
      {0.0086958637, 0.0086958637, 1.0624414, 1.0624414, 2.1043284, 2.1043284},
      {0.010139325, 0.010139325, 1.0625653, 1.0625653, 2.1025118, 2.1025118},
      {0.011580745, 0.011580745, 1.0626571, 1.0626571, 2.1007905, 2.1007905},
      {0.01302009, 0.01302009, 1.0627246, 1.0627246, 2.0991419, 2.0991419},
      {0.014457343, 0.014457343, 1.0627739, 1.0627739, 2.0975492, 2.0975492}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {22.292842, 22.292842, 30.000281, 30.000281, 37.708414, 37.708414},
      {24.060067, 24.060067, 30.000565, 30.000565, 35.941957, 35.941957},
      {25.422176, 25.422176, 30.00082, 30.00082, 34.580297, 34.580297},
      {26.472042, 26.472042, 30.001034, 30.001034, 33.530672, 33.530672},
      {27.281243, 27.281243, 30.001203, 30.001203, 32.721576, 32.721576},
      {27.904947, 27.904947, 30.00133, 30.00133, 32.097889, 32.097889},
      {28.385674, 28.385674, 30.001418, 30.001418, 31.617121, 31.617121},
      {28.756197, 28.756197, 30.001474, 30.001474, 31.246519, 31.246519},
      {29.041778, 29.041778, 30.001502, 30.001502, 30.960835, 30.960835},
      {29.261884, 29.261884, 30.001508, 30.001508, 30.740609, 30.740609}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
