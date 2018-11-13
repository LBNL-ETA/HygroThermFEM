#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;
using MoisThermFEM::State;

class MoistureBC_2D_1 : public testing::Test
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

TEST_F(MoistureBC_2D_1, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    std::vector<double> gridXCoordinates{0,
                                         0.015,
                                         0.025,
                                         0.035,
                                         0.045,
                                         0.055,
                                         0.065,
                                         0.075,
                                         0.085,
                                         0.095,
                                         0.105,
                                         0.115,
                                         0.125,
                                         0.135,
                                         0.15};

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.0;
    const auto initialPressure = 101325.0;

    State state(initialTemperature, initialHumidity, initialPressure, 0);
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
      2050,      /// density
      0.22,      /// porosity
      850,       /// specific heat capacity (dry)
      1.8,       /// thermal conductivity (dry)
      15,        /// diffusion resistance factor (this is mi value)
      {{0, 0},   /// liquid transportation coefficient
                 //{ 27,  1E-8 },
                 //{ 45,  1.1E-8 },
                 //{ 90,  2E-8 },
                 //{ 126, 3.5E-8 },
                 //{ 144, 5E-8 },
                 //{ 162, 1E-7 },
                 //{ 171, 2E-7 },
       {180, 7E-7}},
      {{0, 0},   /// sorption curve
                 // { 0.5,   5.3 },
                 // { 0.65,  8.4 },
                 // { 0.8,   12 },
                 // { 0.93,  17 },
                 // { 0.95,  25 },
                 // { 0.99,  63 },
                 // { 0.995, 83 },
                 // { 0.999, 120 },
       {1, 5.3}});

    MoisThermFEM::MoistureDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (MoisThermFEM::NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    // Create Boundary Conditions
    const auto hc = 20;
    const auto airTemperature = 20;
    const auto airHumidity = 0.5;

    domain.createMoistureBC(1, 2, hc, airHumidity, airTemperature);

    const auto dTime = 36000;
    const auto nSteps = 4;
    /// const auto dTime = 60;
    /// const auto nSteps = 3000;

    auto humidities = NodePool::Instance().properties(MoisThermFEM::Property::humidity);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime);
        auto waterContent = NodePool::Instance().properties(MoisThermFEM::Property::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution{{2.649896,
                                                      0.975243,
                                                      0.483597,
                                                      0.239803,
                                                      0.118912,
                                                      0.058965,
                                                      0.029239,
                                                      0.014499,
                                                      0.007190,
                                                      0.003566,
                                                      0.001770,
                                                      0.000882,
                                                      0.000445,
                                                      0.000236,
                                                      0.000150},
                                                     {2.649960,
                                                      1.410150,
                                                      0.862233,
                                                      0.508375,
                                                      0.292164,
                                                      0.164749,
                                                      0.091549,
                                                      0.050285,
                                                      0.027362,
                                                      0.014777,
                                                      0.007939,
                                                      0.004261,
                                                      0.002315,
                                                      0.001328,
                                                      0.000897},
                                                     {2.649967,
                                                      1.644297,
                                                      1.123855,
                                                      0.737499,
                                                      0.468574,
                                                      0.290061,
                                                      0.175774,
                                                      0.104654,
                                                      0.061399,
                                                      0.035589,
                                                      0.020445,
                                                      0.011710,
                                                      0.006794,
                                                      0.004173,
                                                      0.002975},
                                                     {2.649972,
                                                      1.789028,
                                                      1.307788,
                                                      0.920817,
                                                      0.627800,
                                                      0.416389,
                                                      0.269723,
                                                      0.171207,
                                                      0.106801,
                                                      0.065663,
                                                      0.039940,
                                                      0.024208,
                                                      0.014882,
                                                      0.009701,
                                                      0.007241}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][2 * j], 1e-6);
            // std::cout << solution[i][2 * j] << ",";
        }
        // std::cout << std::endl;
    }
}
