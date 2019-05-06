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
      {0.0014521714, 0.0014521714, 1.060776, 1.060776, 2.1163812, 2.1163812},
      {0.0029038834, 0.0029038834, 1.0613392, 1.0613392, 2.1133332, 2.1133332},
      {0.0043546336, 0.0043546336, 1.0617521, 1.0617521, 2.1106963, 2.1106963},
      {0.0058040514, 0.0058040514, 1.0620571, 1.0620571, 2.1083602, 2.1083602},
      {0.0072518724, 0.0072518724, 1.0622833, 1.0622833, 2.106247, 2.106247},
      {0.0086979139, 0.0086979139, 1.0624517, 1.0624517, 2.1043005, 2.1043005},
      {0.010142054, 0.010142054, 1.0625768, 1.0625768, 2.1024799, 2.1024799},
      {0.011584216, 0.011584216, 1.0626697, 1.0626697, 2.1007549, 2.1007549},
      {0.013024353, 0.013024353, 1.0627382, 1.0627382, 2.0991029, 2.0991029},
      {0.01446244, 0.01446244, 1.0627883, 1.0627883, 2.097507, 2.097507}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
            {22.29179488,22.29179488,30.00106245,30.00106245,37.71890046,37.71890046},
            {24.05870223,24.05870223,30.0026827,30.0026827,35.95862115,35.95862115},
            {25.42105266,25.42105266,30.00443199,30.00443199,34.60036976,34.60036976},
            {26.47154754,26.47154754,30.00609288,30.00609288,33.55240309,33.55240309},
            {27.28161923,27.28161923,30.0075683,30.0075683,32.74388531,32.74388531},
            {27.9063211,27.9063211,30.00882747,30.00882747,32.12013422,32.12013422},
            {28.38808947,28.38808947,30.00987415,30.00987415,31.63894416,31.63894416},
            {28.75963917,28.75963917,30.01072844,30.01072844,31.26774352,31.26774352},
            {29.04619304,29.04619304,30.01141656,30.01141656,30.98139817,30.98139817},
            {29.26719922,29.26719922,30.01196544,30.01196544,30.76051484,30.76051484}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
