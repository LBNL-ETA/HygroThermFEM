#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

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
    const auto liquidPercent = 1.0;

    auto state = HygroThermFEM::State(
      initialTemperature, initialMoistureContent, initialPressure, liquidPercent);
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
            2050,                        /// density
            0.22,                        /// porosity
            850,                         /// specific heat capacity (dry)
            15,                          /// diffusion resistance factor
            {{0.0, 1.8},
             {180, 1.8}},    /// thermal conductivity as function of water content
            {{0,   0},                    /// liquid transportation coefficient
             {27,  1E-8},
             {45,  1.1E-8},
             {90,  2E-8},
             {126, 3.5E-8},
             {144, 5E-8},
             {162, 1E-7},
             {171, 2E-7},
             {180, 7E-7}},
            {{0,     0},   /// sorption curve
             {0.5,   5.3},
             {0.65,  8.4},
             {0.8,   12},
             {0.93,  17},
             {0.95,  25},
             {0.99,  63},
             {0.995, 83},
             {0.999, 120},
             {1,     180}});

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

    /// Create Boundary Conditions
    const auto hc = 1.0;
    const auto airTemperature = 20.0;
    const auto humidity = 0.6;

    domain.createMoistureBCFixedHc(1, 2, airTemperature, hc, humidity);

    const auto dTime = 3600;
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
      {1.089129132, 1.089129132, 0.001302093522, 0.001302093522, 2.981743793e-06, 2.981743793e-06},
      {2.111583133, 2.111583133, 0.003972901399, 0.003972901399, 1.254781245e-05, 1.254781245e-05},
      {2.999134964, 2.999134964, 0.007916291611, 0.007916291611, 3.254787784e-05, 3.254787784e-05},
      {3.72224919, 3.72224919, 0.01293991832, 0.01293991832, 6.675532267e-05, 6.675532267e-05},
      {4.333515338, 4.333515338, 0.01894277529, 0.01894277529, 0.0001190156958, 0.0001190156958},
      {4.861925859, 4.861925859, 0.02585262679, 0.02585262679, 0.0001932689657, 0.0001932689657},
      {5.349567861, 5.349567861, 0.0336129523, 0.0336129523, 0.0002935482315, 0.0002935482315},
      {6.150137728, 6.150137728, 0.04217641849, 0.04217641849, 0.0004239662948, 0.0004239662948},
      {6.864459538, 6.864459538, 0.05150164249, 0.05150164249, 0.0005886945982, 0.0005886945982},
      {7.511602342, 7.511602342, 0.06157470173, 0.06157470173, 0.000792020886, 0.000792020886}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());


    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.315929677, 1.315929677, 0.6840152658, 0.6840152658, 0.5119375064, 0.5119375064},
      {2.163637359, 2.163637359, 1.430334777, 1.430334777, 1.199293069, 1.199293069},
      {2.877362614, 2.877362614, 2.151732423, 2.151732423, 1.912125752, 1.912125752},
      {3.533487377, 3.533487377, 2.837848324, 2.837848324, 2.604961318, 2.604961318},
      {4.152881584, 4.152881584, 3.489969466, 3.489969466, 3.267323099, 3.267323099},
      {4.7414041, 4.7414041, 4.110206843, 4.110206843, 3.898155439, 3.898155439},
      {5.301411772, 5.301411772, 4.700377755, 4.700377755, 4.498552903, 4.498552903},
      {5.834359222, 5.834359222, 5.262010079, 5.262010079, 5.069934374, 5.069934374},
      {6.341119777, 6.341119777, 5.796248026, 5.796248026, 5.613513172, 5.613513172},
      {6.823006749, 6.823006749, 6.304365864, 6.304365864, 6.130548244, 6.130548244}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
