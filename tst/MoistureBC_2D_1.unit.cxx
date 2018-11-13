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
        auto & node1 = MoisThermFEM::NodePool::Instance().getNode(2 * i + 1);
        auto & node2 = MoisThermFEM::NodePool::Instance().getNode(2 * i + 2);
        auto & node3 = MoisThermFEM::NodePool::Instance().getNode(2 * i);
        auto & node4 = MoisThermFEM::NodePool::Instance().getNode(2 * i - 1);
        domain.createElement(node1, node2, node3, node4, material);
    }

    // Create Boundary Conditions
    const auto hc = 20;
    const auto airTemperature = 20;
    const auto airHumidity = 0.5;

    auto & node1 = MoisThermFEM::NodePool::Instance().getNode(1);
    auto & node2 = MoisThermFEM::NodePool::Instance().getNode(2);

    domain.createMoistureBC(node1, node2, hc, airHumidity, airTemperature);

    const auto dTime = 36000;
    const auto nSteps = 4;
    /// const auto dTime = 60;
    /// const auto nSteps = 3000;

    auto humidities = NodePool::Instance().nodeProperties(MoisThermFEM::StateProperty::humidity);
    std::vector<std::vector<double>> solution;

    for(unsigned i = 0; i < nSteps; ++i)
    {
        humidities = domain.transient(humidities, dTime);
        solution.push_back(material.waterContent(humidities));
    }

    std::vector<std::vector<double>> correctSolution{{2.6498955,
                                                      0.97524335,
                                                      0.48359673,
                                                      0.23980252,
                                                      0.11891159,
                                                      0.058965089,
                                                      0.029239302,
                                                      0.014499207,
                                                      0.0071902274,
                                                      0.0035663744,
                                                      0.0017703551,
                                                      0.00088167597,
                                                      0.00044487237,
                                                      0.00023607421,
                                                      0.00014973783},
                                                     {2.6499595,
                                                      1.4101503,
                                                      0.86223326,
                                                      0.50837455,
                                                      0.29216405,
                                                      0.16474869,
                                                      0.091549415,
                                                      0.050285239,
                                                      0.027362077,
                                                      0.014777367,
                                                      0.0079385003,
                                                      0.0042609241,
                                                      0.0023152761,
                                                      0.001328245,
                                                      0.00089724477},
                                                     {2.6499672,
                                                      1.6442965,
                                                      1.1238547,
                                                      0.73749888,
                                                      0.46857355,
                                                      0.29006141,
                                                      0.17577438,
                                                      0.10465423,
                                                      0.061399201,
                                                      0.03558884,
                                                      0.020444751,
                                                      0.011710349,
                                                      0.0067939179,
                                                      0.0041728745,
                                                      0.0029749206},
                                                     {2.6499719,
                                                      1.7890282,
                                                      1.3077879,
                                                      0.92081676,
                                                      0.6277995,
                                                      0.41638853,
                                                      0.26972258,
                                                      0.17120685,
                                                      0.10680059,
                                                      0.065663392,
                                                      0.039939966,
                                                      0.024208203,
                                                      0.014881826,
                                                      0.0097006528,
                                                      0.0072409374}};

    EXPECT_EQ(solution.size(), correctSolution.size());

    for(auto i = 0u; i < correctSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctSolution[i][j], solution[i][2 * j], 1e-6);
        }
    }
}
