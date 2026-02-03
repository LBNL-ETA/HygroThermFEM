#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::SimulationProperties;

class MultiDomain_2D_ThermalConductivityDependent_1 : public testing::Test
{
protected:
    void SetUp() override
    {
        const auto relaxationParameter{0.8};
        const auto errorTolerance{1e-5};
        const auto numberOfIterations{20u};
        SimulationProperties::Instance().setIterationParameters(
          relaxationParameter, errorTolerance, numberOfIterations);
    }

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
        SimulationProperties::Instance().reset();
    }
};

TEST_F(MultiDomain_2D_ThermalConductivityDependent_1, TestExample_1)
{
    const auto excludeWaterLiquidTransportation{false};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{false};
    const auto thermalConductivityMoistureAndTemperatureDependent{true};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    std::vector<double> gridXCoordinates{0, 0.05, 0.1};

    const double initialTemperature = 0.0;
    const double initialMoistureContent = 0.99;
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

    // Material Properties (Cottaer Sandstone - Thermal Conductivity Dependent on Temperature and
    // Moisture)
    constexpr double thermalConductivityDry{1.8};
    constexpr double density{2050.0};
    constexpr double porosity{0.22};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 2.5}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 3.1}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
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
    const auto hc = 5.0;
    const auto airTemperature = 10.0;
    const auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    domain.createBC_FixedHc(1, 2, bcCoeff);
    domain.createBC_FixedHc(5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 24;

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
        {13.313624,13.313624,25.249588,25.249588,13.313624,13.313624},
        {12.614011,12.614011,25.243483,25.243483,12.614011,12.614011},
        {11.944162,11.944162,25.236386,25.236386,11.944162,11.944162},
        {11.503589,11.503589,25.228285,25.228285,11.503589,11.503589},
        {11.062118,11.062118,25.219166,25.219166,11.062118,11.062118},
        {10.620338,10.620338,25.209013,25.209013,10.620338,10.620338},
        {10.178964,10.178964,25.197814,25.197814,10.178964,10.178964},
        {9.742766,9.742766,25.185568,25.185568,9.742766,9.742766},
        {9.312022,9.312022,25.172271,25.172271,9.312022,9.312022},
        {8.886953,8.886953,25.157921,25.157921,8.886953,8.886953},
        {8.467741,8.467741,25.142513,25.142513,8.467741,8.467741},
        {8.102511,8.102511,25.126044,25.126044,8.102511,8.102511},
        {7.751963,7.751963,25.10851,25.10851,7.751963,7.751963},
        {7.406777,7.406777,25.089907,25.089907,7.406777,7.406777},
        {7.067021,7.067021,25.070232,25.070232,7.067021,7.067021},
        {6.732639,6.732639,25.049481,25.049481,6.732639,6.732639},
        {6.403474,6.403474,25.027649,25.027649,6.403474,6.403474},
        {6.079589,6.079589,25.004732,25.004732,6.079589,6.079589},
        {5.76104,5.76104,24.991884,24.991884,5.76104,5.76104},
        {5.447875,5.447875,24.981316,24.981316,5.447875,5.447875},
        {5.218007,5.218007,24.970287,24.970287,5.218007,5.218007},
        {5.062973,5.062973,24.958795,24.958795,5.062973,5.062973},
        {4.910761,4.910761,24.946839,24.946839,4.910761,4.910761},
        {4.761387,4.761387,24.934418,24.934418,4.761387,4.761387}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
        {0.346934,0.346934,0.269692,0.269692,0.346934,0.346934},
        {0.607114,0.607114,0.542994,0.542994,0.607114,0.607114},
        {0.835532,0.835532,0.78634,0.78634,0.835532,0.835532},
        {1.042068,1.042068,1.003079,1.003079,1.042068,1.042068},
        {1.231444,1.231444,1.199349,1.199349,1.231444,1.231444},
        {1.407421,1.407421,1.380147,1.380147,1.407421,1.407421},
        {1.573099,1.573099,1.549307,1.549307,1.573099,1.573099},
        {1.730687,1.730687,1.70951,1.70951,1.730687,1.730687},
        {1.881956,1.881956,1.862797,1.862797,1.881956,1.881956},
        {2.028299,2.028299,2.01073,2.01073,2.028299,2.028299},
        {2.170795,2.170795,2.154507,2.154507,2.170795,2.170795},
        {2.310222,2.310222,2.295057,2.295057,2.310222,2.310222},
        {2.447309,2.447309,2.433027,2.433027,2.447309,2.447309},
        {2.582514,2.582514,2.568981,2.568981,2.582514,2.582514},
        {2.716212,2.716212,2.703326,2.703326,2.716212,2.716212},
        {2.848693,2.848693,2.836374,2.836374,2.848693,2.848693},
        {2.980192,2.980192,2.968375,2.968375,2.980192,2.980192},
        {3.11088,3.11088,3.099513,3.099513,3.11088,3.11088},
        {3.240877,3.240877,3.229917,3.229917,3.240877,3.240877},
        {3.370266,3.370266,3.359679,3.359679,3.370266,3.370266},
        {3.499082,3.499082,3.488841,3.488841,3.499082,3.499082},
        {3.627343,3.627343,3.617423,3.617423,3.627343,3.627343},
        {3.75507,3.75507,3.745451,3.745451,3.75507,3.75507},
        {3.882267,3.882267,3.87293,3.87293,3.882267,3.882267}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
