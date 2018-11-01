#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class MoistureBC_2D_2 : public testing::Test
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

TEST_F(MoistureBC_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    std::vector<double> gridXCoordinates{0.15, 0.05, 0.00};

    MoisThermFEM::State state(293.15, 0, 101325, 0);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
    }

    auto & material = MaterialPool::Instance().createMaterial(
      "Cottaer Sandstone",
      2050,      /// density
      0.22,      /// porosity
      850,       /// specific heat capacity (dry)
      1.8,       /// thermal conductivity (dry)
      15,        /// diffusion resistance factor (this is mi value)
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

    MoisThermFEM::MoistureDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        auto node1 = NodePool::Instance().getNode(2 * i + 1);
        auto node2 = NodePool::Instance().getNode(2 * i + 2);
        auto node3 = NodePool::Instance().getNode(2 * i);
        auto node4 = NodePool::Instance().getNode(2 * i - 1);
        domain.createElement(node2, node3, node4, node1, material);
    }

    // Create Boundary Conditions
    const auto hc = 20;
    const auto airTemperature = 293.15;
    const auto humidity = 0.5;

    auto node1 = NodePool::Instance().getNode(5);
    auto node2 = NodePool::Instance().getNode(6);

    domain.createMoistureBC(node1, node2, hc, humidity, airTemperature);

    const auto dTime = 36000;
    const auto nSteps = 4;

    auto humidities = NodePool::Instance().nodeProperties(MoisThermFEM::Property::humidity);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime);
        solution.push_back(material.waterContent(humidities));
    }

    std::vector<std::vector<double>> correctSolution{
      {0.0025399387, 0.0025399387, 0.13271657, 0.13271657, 5.2992847, 5.2992847},
      {0.0074773153, 0.0074773153, 0.26052713, 0.26052713, 5.2999506, 5.2999506},
      {0.014676106, 0.014676106, 0.38362762, 0.38362762, 5.2999519, 5.2999519},
      {0.024006785, 0.024006785, 0.50222163, 0.50222163, 5.299953, 5.299953}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
        }
    }
}
