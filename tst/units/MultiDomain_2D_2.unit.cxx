#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;

class MultiDomain_2D_2 : public testing::Test
{
public:
    HygroThermFEM::MultiDomain domain;

    const double dTime{3600};
    const size_t nSteps{10u};

    const double initialTemperature{0.0};
    const double initialMoistureContent{0.0};
    const double initialPressure{101325};
    const double liquidPercent{1.0};

    const double hc{1.0};
    const double airTemperature{20.0};
    const double humidity{0.6};

protected:
    void SetUp() override
    {
        std::vector<double> gridXCoordinates{0, 0.05, 0.1};

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

        // Material Properties (Cottaer Sandstone)
        const double thermalConductivityDry{1.8};
        const double density{2050.0};
        const double porosity{0.22};
        const double specificHeatCapacityDry{850.0};
        const double diffusionResistanceFactor{15.0};
        const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
          {0.0, 1.8}, {180, 1.8}};
        const double thermalConductivityMeasuredAtTemperature{0};
        const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
          {0.0, 1.8}, {1, 1.8}};
        const double thermalConductivityMeasuredAtHumidity{0};
        const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                                  {27, 1E-8},
                                                                                  {45, 1.1E-8},
                                                                                  {90, 2E-8},
                                                                                  {126, 3.5E-8},
                                                                                  {144, 5E-8},
                                                                                  {162, 1E-7},
                                                                                  {171, 2E-7},
                                                                                  {180, 7E-7}};

        const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                                {0.5, 5.3},
                                                                                {0.65, 8.4},
                                                                                {0.8, 12},
                                                                                {0.93, 17},
                                                                                {0.95, 25},
                                                                                {0.99, 63},
                                                                                {0.995, 83},
                                                                                {0.999, 120},
                                                                                {1, 180}};

        auto & material =
          MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone",
                                                       thermalConductivityDry,
                                                       density,
                                                       porosity,
                                                       specificHeatCapacityDry,
                                                       diffusionResistanceFactor,
                                                       thermalConductivityMoistureDependent,
                                                       thermalConductivityMeasuredAtTemperature,
                                                       thermalConductivityTemperatureDependent,
                                                       thermalConductivityMeasuredAtHumidity,
                                                       liquidTransportationCurve,
                                                       moistureStorageFunction);

        /// Create elements
        for(size_t i = 1; i <= (HygroThermFEM::maxNodeIndex() - 2) / 2; ++i)
        {
            const auto node1 = 2u * i + 1u;
            const auto node2 = 2u * i + 2u;
            const auto node3 = 2u * i;
            const auto node4 = 2u * i - 1u;
            createElement(domain, node1, node2, node3, node4, material.name());
        }

        const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

        createBC_FixedHc(domain, 1, 2, bcCoeff);
    }

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
    }
};

TEST_F(MultiDomain_2D_2, Substitution)
{
    SCOPED_TRACE("Begin Test: Simple two elements example with moisture and heat transfer.");

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    HygroThermFEM::TransientSubstitutionSolver solver{domain};
    for(auto i = 0u; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {1.089281, 1.089281, 0.001246, 0.001246, 0.000003, 0.000003},
      {2.112231, 2.112231, 0.003780, 0.003780, 0.000012, 0.000012},
      {3.000407, 3.000407, 0.007522, 0.007522, 0.000030, 0.000030},
      {3.724351, 3.724351, 0.012294, 0.012294, 0.000062, 0.000062},
      {4.336618, 4.336618, 0.018008, 0.018008, 0.000111, 0.000111},
      {4.866186, 4.866186, 0.024596, 0.024596, 0.000180, 0.000180},
      {5.360410, 5.360410, 0.032009, 0.032009, 0.000274, 0.000274},
      {6.163777, 6.163777, 0.040203, 0.040203, 0.000396, 0.000396},
      {6.881329, 6.881329, 0.049142, 0.049142, 0.000550, 0.000550},
      {7.531928, 7.531928, 0.058814, 0.058814, 0.000740, 0.000740},
    };

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.278962, 1.278962, 0.664542, 0.664542, 0.497362, 0.497362},
      {2.103073, 2.103073, 1.389887, 1.389887, 1.165354, 1.165354},
      {2.799477, 2.799477, 2.092455, 2.092455, 1.859221, 1.859221},
      {3.443158, 3.443158, 2.763086, 2.763086, 2.535696, 2.535696},
      {4.052248, 4.052248, 3.402329, 3.402329, 3.184303, 3.184303},
      {4.631946, 4.631946, 4.011701, 4.011701, 3.803543, 3.803543},
      {5.184334, 5.184334, 4.592612, 4.592612, 4.394094, 4.394094},
      {5.710388, 5.710388, 5.146121, 5.146121, 4.956919, 4.956919},
      {6.211406, 6.211406, 5.673433, 5.673433, 5.493161, 5.493161},
      {6.688342, 6.688342, 6.175609, 6.175609, 6.003903, 6.003903},
    };

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}

TEST_F(MultiDomain_2D_2, Substitution_Repeat)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    HygroThermFEM::TransientSubstitutionSolver solver{domain};
    for(auto i = 0u; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {1.089281, 1.089281, 0.001246, 0.001246, 0.000003, 0.000003},
      {2.112231, 2.112231, 0.003780, 0.003780, 0.000012, 0.000012},
      {3.000407, 3.000407, 0.007522, 0.007522, 0.000030, 0.000030},
      {3.724351, 3.724351, 0.012294, 0.012294, 0.000062, 0.000062},
      {4.336618, 4.336618, 0.018008, 0.018008, 0.000111, 0.000111},
      {4.866186, 4.866186, 0.024596, 0.024596, 0.000180, 0.000180},
      {5.360410, 5.360410, 0.032009, 0.032009, 0.000274, 0.000274},
      {6.163777, 6.163777, 0.040203, 0.040203, 0.000396, 0.000396},
      {6.881329, 6.881329, 0.049142, 0.049142, 0.000550, 0.000550},
      {7.531928, 7.531928, 0.058814, 0.058814, 0.000740, 0.000740},
    };

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.278962, 1.278962, 0.664542, 0.664542, 0.497362, 0.497362},
      {2.103073, 2.103073, 1.389887, 1.389887, 1.165354, 1.165354},
      {2.799477, 2.799477, 2.092455, 2.092455, 1.859221, 1.859221},
      {3.443158, 3.443158, 2.763086, 2.763086, 2.535696, 2.535696},
      {4.052248, 4.052248, 3.402329, 3.402329, 3.184303, 3.184303},
      {4.631946, 4.631946, 4.011701, 4.011701, 3.803543, 3.803543},
      {5.184334, 5.184334, 4.592612, 4.592612, 4.394094, 4.394094},
      {5.710388, 5.710388, 5.146121, 5.146121, 4.956919, 4.956919},
      {6.211406, 6.211406, 5.673433, 5.673433, 5.493161, 5.493161},
      {6.688342, 6.688342, 6.175609, 6.175609, 6.003903, 6.003903},
    };

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}

TEST_F(MultiDomain_2D_2, Sundials)
{
    SCOPED_TRACE("Begin Test: Repeatability test.");

    auto temperatures{properties(HygroThermFEM::Variable::temperature)};
    auto humidities{properties(HygroThermFEM::Variable::humidity)};
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    Sundials::SolverIDA solver{domain};
    for(auto i = 0u; i < nSteps; ++i)
    {
        auto aSolution = solver.transient(temperatures, humidities, dTime, i);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {1.133346, 1.133346, 0.000632, 0.000632, 0.000000, 0.000000},
      {2.270425, 2.270425, 0.002742, 0.002742, 0.000004, 0.000004},
      {3.319572, 3.319572, 0.006346, 0.006346, 0.000016, 0.000016},
      {4.015458, 4.015458, 0.010335, 0.010335, 0.000034, 0.000034},
      {4.699985, 4.699985, 0.016157, 0.016157, 0.000069, 0.000069},
      {5.210962, 5.210962, 0.022216, 0.022216, 0.000118, 0.000118},
      {6.066415, 6.066415, 0.029659, 0.029659, 0.000192, 0.000192},
      {6.855831, 6.855831, 0.037585, 0.037585, 0.000288, 0.000288},
      {7.601706, 7.601706, 0.046844, 0.046844, 0.000421, 0.000421},
      {8.238425, 8.238425, 0.056369, 0.056369, 0.000580, 0.000580}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < correctWaterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctWaterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {1.421933, 1.421933, 0.661658, 0.661658, 0.423883, 0.423883},
      {2.153662, 2.153662, 1.408789, 1.408789, 1.160040, 1.160040},
      {2.820270, 2.820270, 2.113518, 2.113518, 1.876861, 1.876861},
      {3.463315, 3.463315, 2.789386, 2.789386, 2.563932, 2.563932},
      {4.070768, 4.070768, 3.429730, 3.429730, 3.215379, 3.215379},
      {4.652676, 4.652676, 4.041034, 4.041034, 3.836678, 3.836678},
      {5.203143, 5.203143, 4.620478, 4.620478, 4.425917, 4.425917},
      {5.735399, 5.735399, 5.180103, 5.180103, 4.994835, 4.994835},
      {6.235458, 6.235458, 5.706820, 5.706820, 5.530549, 5.530549},
      {6.709401, 6.709401, 6.205772, 6.205772, 6.037964, 6.037964},
    };

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < correctTemperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < correctTemperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}