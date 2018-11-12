#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;
using MoisThermFEM::State;

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

    auto & material =
      MaterialPool::Instance().createMaterial("Cottaer Sandstone",
                                              2050,      /// density
                                              0.22,      /// porosity
                                              850,       /// specific heat capacity (dry)
                                              1.8,       /// thermal conductivity (dry)
                                              15,        /// diffusion resistance factor
                                              {{0, 0},   /// liquid transportation coefficient
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
        auto & node1 = NodePool::Instance().getNode(2 * i + 1);
        auto & node2 = NodePool::Instance().getNode(2 * i + 2);
        auto & node3 = NodePool::Instance().getNode(2 * i);
        auto & node4 = NodePool::Instance().getNode(2 * i - 1);
        domain.createElement(node1, node2, node3, node4, material);
    }

    const auto dTime = 360;
    /// const auto dTime = 3600;
    const auto nSteps = 10;
    /// const auto nSteps = 10;

    auto temperatures = NodePool::Instance().nodeProperties(MoisThermFEM::Property::temperature);
    auto humidities = NodePool::Instance().nodeProperties(MoisThermFEM::Property::humidity);
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
      {0.0020830641, 0.0020830641, 1.0610351, 1.0610351, 2.1150155, 2.1150155},
      {0.0041649337, 0.0041649337, 1.0617831, 1.0617831, 2.110805, 2.110805},
      {0.0062445946, 0.0062445946, 1.0623278, 1.0623278, 2.1071533, 2.1071533},
      {0.0083213149, 0.0083213149, 1.0627264, 1.0627264, 2.103911, 2.103911},
      {0.010394582, 0.010394582, 1.0630185, 1.0630185, 2.1009727, 2.1009727},
      {0.012464051, 0.012464051, 1.0632321, 1.0632321, 2.0982621, 2.0982621},
      {0.014529497, 0.014529497, 1.0633874, 1.0633874, 2.0957236, 2.0957236},
      {0.016590782, 0.016590782, 1.063499, 1.063499, 2.0933161, 2.0933161},
      {0.018647828, 0.018647828, 1.0635777, 1.0635777, 2.0910089, 2.0910089},
      {0.0207006, 0.0207006, 1.0636316, 1.0636316, 2.088779, 2.088779}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    std::cout << "Water Content" << std::endl;

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            // EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
            std::cout << waterContentSolution[i][j] << ",";
        }
        std::cout << std::endl;
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {302.3175172, 302.3175172, 310.0016490, 310.0016490, 317.6932556, 317.6932556},
      {304.1027373, 304.1027373, 310.0035788, 310.0035788, 315.9149897, 315.9149897},
      {305.4780085, 305.4780085, 310.0054784, 310.0054784, 314.5442326, 314.5442326},
      {306.5375218, 306.5375218, 310.0072002, 310.0072002, 313.4876682, 313.4876682},
      {307.3538038, 307.3538038, 310.0086886, 310.0086886, 312.6733264, 312.6733264},
      {307.9827088, 307.9827088, 310.0099367, 310.0099367, 312.0457085, 312.0457085},
      {308.4672557, 308.4672557, 310.0109618, 310.0109618, 311.5620227, 311.5620227},
      {308.8405817, 308.8405817, 310.0117916, 310.0117916, 311.1892785, 311.1892785},
      {309.1282137, 309.1282137, 310.0124561, 310.0124561, 310.9020434, 310.9020434},
      {309.3498179, 309.3498179, 310.0129841, 310.0129841, 310.6807131, 310.6807131}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    std::cout << "Temperatures" << std::endl;

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            // EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
            std::cout << temperatureSolution[i][j] << ",";
        }
        std::cout << std::endl;
    }
}
