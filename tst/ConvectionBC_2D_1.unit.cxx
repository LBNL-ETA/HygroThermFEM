#include <gtest/gtest.h>
#include <memory>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class ConvectionBC_2D_1 : public testing::Test
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

TEST_F(ConvectionBC_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with simple conduction.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    auto & nodePool = NodePool::Instance();
    auto & materialPool = MaterialPool::Instance();

    const auto node1 = nodePool.createNode(1, 15, 5);
    const auto node2 = nodePool.createNode(2, 15, 0);
    const auto node3 = nodePool.createNode(3, 5, 5);
    const auto node4 = nodePool.createNode(4, 5, 0);
    const auto node5 = nodePool.createNode(5, 0, 5);
    const auto node6 = nodePool.createNode(6, 0, 0);

    auto & material = materialPool.createMaterial("Test Material",
                                                  2050,      /// Density
                                                  0.00,      /// Porosity
                                                  850,       /// Specific Heat Capacity (dry)
                                                  1,         /// Thermal Conductivity (dry)
                                                  15E-6,     /// Diffusion Resistance Factor
                                                  {{0, 0},   /// Liquid Transportation Coefficient
                                                   {27, 1E-8},
                                                   {45, 1.1E-8},
                                                   {90, 2E-8},
                                                   {126, 3.5E-8},
                                                   {144, 5E-8},
                                                   {162, 1E-7},
                                                   {171, 2E-7},
                                                   {180, 7E-7}},
                                                  {{0, 0},   /// Moisture Storage Function
                                                   {0.5, 5.3},
                                                   {0.65, 8.4},
                                                   {0.8, 12},
                                                   {0.93, 17},
                                                   {0.95, 25},
                                                   {0.99, 63},
                                                   {0.995, 83},
                                                   {0.999, 120},
                                                   {1, 180}});

    MoisThermFEM::ThermalDomain domain;

    domain.createElement(node3, node4, node2, node1, material);
    domain.createElement(node6, node4, node3, node5, material);

    // Create Boundary Conditions
    const auto hc1 = 20.0;
    const auto temperatureAir1 = -18.0;

    const auto hc2 = 2.4;
    const auto temperatureAir2 = 21.0;

    domain.createConvectionBC(node1, node2, hc1, temperatureAir1);
    domain.createConvectionBC(node6, node5, hc2, temperatureAir2);

    auto solution = domain.steadyState();

    std::vector<double> correctSolution = {
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        EXPECT_NEAR(correctSolution[i], solution[i], 1e-6);
    }
}
