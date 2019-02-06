#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MoistureBC_2D_3 : public testing::Test
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

TEST_F(MoistureBC_2D_3, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture transfer.");

    std::vector<double> gridXCoordinates{0, 0.05, 0.1, 0.15};

    const double initialTemperature = 20;
    const double initialMoistureContent = 0;
    const double initialPressure = 0;

    HygroThermFEM::State state(initialTemperature, initialMoistureContent, initialPressure, 0);
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
      {{0.0, 1.8}, {180, 1.8}},   /// thermal conductivity (dry)
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

    HygroThermFEM::MoistureDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    // Create Boundary Conditions
    const auto airTemperature = 20;
    const auto humidity = 0.5;

    domain.createMoistureBCVariableHc(1, 2, humidity, airTemperature);

    const auto dTime = 36000;
    const auto nSteps = 4;

    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime).solution;
        auto waterContent = NodePool::Instance().properties(HygroThermFEM::Variable::water);
        solution.push_back(waterContent);
    }

    std::vector<std::vector<double>> correctSolution = {{4.56238485,
                                                         4.56238485,
                                                         0.165365838,
                                                         0.165365838,
                                                         0.00600161797,
                                                         0.00600161797,
                                                         0.000434491517,
                                                         0.000434491517},
                                                        {5.15447905,
                                                         5.15447905,
                                                         0.340625261,
                                                         0.340625261,
                                                         0.0179586848,
                                                         0.0179586848,
                                                         0.00170316829,
                                                         0.00170316829},
                                                        {5.23276638,
                                                         5.23276638,
                                                         0.506652687,
                                                         0.506652687,
                                                         0.0351475929,
                                                         0.0351475929,
                                                         0.00412440184,
                                                         0.00412440184},
                                                        {5.24446337,
                                                         5.24446337,
                                                         0.661869378,
                                                         0.661869378,
                                                         0.0568487783,
                                                         0.0568487783,
                                                         0.00794142158,
                                                         0.00794142158}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    std::cout.precision(10);
    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        std::cout << "{";
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            //EXPECT_NEAR(correctSolution[i][j], solution[i][j], 1e-6);
            std::cout << solution[i][j] << ",";
        }
        std::cout << "}," << std::endl;
    }
}
