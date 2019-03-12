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
                                                              2050, /// density
                                                              0.22, /// porosity
                                                              850,  /// specific heat capacity (dry)
                                                              15,   /// diffusion resistance factor
                                                              {{0.0, 1.8},
                                                              {180, 1.8}}, /// thermal conductivity as function of water content
                                                              {{0,   0},   /// liquid transportation coefficient
                                                              {27,  1E-8},
                                                              {45,  1.1E-8},
                                                              {90,  2E-8},
                                                              {126, 3.5E-8},
                                                              {144, 5E-8},
                                                              {162, 1E-7},
                                                              {171, 2E-7},
                                                              {180, 7E-7}},
                                                              {{0,     0}, /// sorption curve
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
      {0.001452171568, 0.001452171568, 1.060775964, 1.060775964, 2.116381217, 2.116381217},
      {0.002903884058, 0.002903884058, 1.061339178, 1.061339178, 2.113333192, 2.113333192},
      {0.004354634927, 0.004354634927, 1.061752111, 1.061752111, 2.110696288, 2.110696288},
      {0.005804053712, 0.005804053712, 1.06205707, 1.06205707, 2.108360188, 2.108360188},
      {0.007251875967, 0.007251875967, 1.062283359, 1.062283359, 2.106246946, 2.106246946},
      {0.008697918973, 0.008697918973, 1.062451669, 1.062451669, 2.104300483, 2.104300483},
      {0.01014206097, 0.01014206097, 1.062576855, 1.062576855, 2.102479858, 2.102479858},
      {0.01158422445, 0.01158422445, 1.062669726, 1.062669726, 2.100754828, 2.100754828},
      {0.01302436316, 0.01302436316, 1.06273824, 1.06273824, 2.099102843, 2.099102843},
      {0.01446245241, 0.01446245241, 1.06278831, 1.06278831, 2.097506969, 2.097506969}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-8);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {22.29180285, 22.29180285, 30.00106456, 30.00106456, 37.71891899, 37.71891899},
      {24.05871516, 24.05871516, 30.0026883, 30.0026883, 35.95864874, 35.95864874},
      {25.42106882, 25.42106882, 30.00444131, 30.00444131, 34.60040144, 34.60040144},
      {26.47156588, 26.47156588, 30.00610562, 30.00610562, 33.55243628, 33.55243628},
      {27.28163911, 27.28163911, 30.00758398, 30.00758398, 32.74391875, 32.74391875},
      {27.9063421, 27.9063421, 30.00884556, 30.00884556, 32.12016732, 32.12016732},
      {28.38811132, 28.38811132, 30.00989419, 30.00989419, 31.6389767, 31.6389767},
      {28.75966167, 28.75966167, 30.01075002, 30.01075002, 31.26777545, 31.26777545},
      {29.04621604, 29.04621604, 30.01143935, 30.01143935, 30.98142954, 30.98142954},
      {29.26722262, 29.26722262, 30.01198917, 30.01198917, 30.7605457, 30.7605457}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
