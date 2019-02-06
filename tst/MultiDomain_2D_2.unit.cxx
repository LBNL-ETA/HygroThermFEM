#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using MoisThermFEM::NodePool;
using MoisThermFEM::MaterialPool;

class MultiDomain_2D_2 : public testing::Test
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

TEST_F(MultiDomain_2D_2, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;

    auto state =
      MoisThermFEM::State(initialTemperature, initialMoistureContent, initialPressure, 0);
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

    MoisThermFEM::MultiDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (NodePool::Instance().maxIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        domain.createElement(node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    const auto hc = 1.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.2;

    domain.createMoistureBCFixedHc(1, 2, airTemperature, hc, humidity);

    const auto dTime = 3600;
    const auto nSteps = 10;

    auto temperatures = NodePool::Instance().properties(MoisThermFEM::Variable::temperature);
    auto humidities = NodePool::Instance().properties(MoisThermFEM::Variable::humidity);
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
        {0.4420797505,0.4420797505,0.000505285817,0.000505285817,1.130261458e-006,1.130261458e-006},
        {0.856271176,0.856271176,0.001512934644,0.001512934644,4.602065313e-006,4.602065313e-006},
        {1.242976504,1.242976504,0.003013152236,0.003013152236,1.169403466e-005,1.169403466e-005},
        {1.602696465,1.602696465,0.004994365415,0.004994365415,2.374398005e-005,2.374398005e-005},
        {1.93601288,1.93601288,0.007443537713,0.007443537713,4.214214088e-005,4.214214088e-005},
        {2.243590435,2.243590435,0.01034616053,0.01034616053,6.832391812e-005,6.832391812e-005},
        {2.526173038,2.526173038,0.01368633634,0.01368633634,0.0001037621799,0.0001037621799},
        {2.784576474,2.784576474,0.01744692422,0.01744692422,0.0001499593752,0.0001499593752},
        {3.019679526,3.019679526,0.02160971038,0.02160971038,0.0002084396125,0.0002084396125},
        {3.232414572,3.232414572,0.02615558638,0.02615558638,0.0002807408105,0.0002807408105}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
        {0.6580943586,0.6580943586,0.342072536,0.342072536,0.2560173255,0.2560173255},
        {1.101397422,1.101397422,0.7255088229,0.7255088229,0.6073988074,0.6073988074},
        {1.4910601,1.4910601,1.107988575,1.107988575,0.9820551517,0.9820551517},
        {1.862567238,1.862567238,1.483741583,1.483741583,1.357532254,1.357532254},
        {2.224183592,2.224183592,1.852111335,1.852111335,1.727689981,1.727689981},
        {2.577977705,2.577977705,2.213140777,2.213140777,2.091015802,2.091015802},
        {2.924556759,2.924556759,2.566962748,2.566962748,2.447228631,2.447228631},
        {3.264173984,3.264173984,2.913718516,2.913718516,2.796363484,2.796363484},
        {3.596994241,3.596994241,3.253547661,3.253547661,3.138533635,3.138533635},
        {3.923158769,3.923158769,3.5865871,3.5865871,3.473870057,3.473870057}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
