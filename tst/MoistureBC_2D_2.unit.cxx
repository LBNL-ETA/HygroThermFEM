#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

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

    const auto domainTemperature = 20.0;
    const auto domainHumidity = 0.0;
    const auto domainPressure = 101325;
    const auto liquidPercentage = 1.0;

    const HygroThermFEM::State state(
      domainTemperature, domainHumidity, domainPressure, liquidPercentage);
    size_t nodeIndex = 0;
    for(auto val : gridXCoordinates)
    {
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.05, state);
        ++nodeIndex;
        NodePool::Instance().createNode(nodeIndex, val, 0.00, state);
    }

    auto & material = MaterialPool::Instance().createSolidMaterial(
      "Cottaer Sandstone",
      2050,                         /// density
      0.22,                         /// porosity
      850,                          /// specific heat capacity (dry)
      15,                           /// diffusion resistance factor (this is mi value)
      {{0.0, 1.8}, {180.0, 1.8}},   /// thermal conductivity (dry)
      {{0, 0},                      /// liquid transportation coefficient
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

    HygroThermFEM::MoistureDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node2, node3, node4, node1, material.name());
    }

    // Create Boundary Conditions
    const auto ambientTemperature = 20;
    const auto ambientHumidity = 0.2;

    domain.createMoistureBCVariableHc(5, 6, ambientHumidity, ambientTemperature);

    const auto dTime = 3600;
    const auto nSteps = 4;

    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> solution;

    for(size_t i = 0u; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime).solution;
        auto waterContent = NodePool::Instance().properties(HygroThermFEM::Variable::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution{
      {3.748586e-06, 3.748586e-06, 0.001924969331, 0.001924969331, 0.7428227662, 0.7428227662},
      {1.364027149e-05, 1.364027149e-05, 0.005083314705, 0.005083314705, 1.221651619, 1.221651619},
      {3.11968812e-05, 3.11968812e-05, 0.009029288908, 0.009029288908, 1.530317721, 1.530317721},
      {5.737783803e-05, 5.737783803e-05, 0.01347560965, 0.01347560965, 1.729301858, 1.729301858}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-9);
        }
    }
}
