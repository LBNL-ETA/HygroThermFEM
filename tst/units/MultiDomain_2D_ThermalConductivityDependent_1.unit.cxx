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
    const double thermalConductivityDry{1.8};
    const double density{2050.0};
    const double porosity{0.22};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 2.5}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 3.1}};
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

    HygroThermFEM::MultiDomain domain;

    /// Create elements
    for(size_t i = 1; i <= (HygroThermFEM::maxNodeIndex() - 2) / 2; ++i)
    {
        const auto node1 = 2u * i + 1u;
        const auto node2 = 2u * i + 2u;
        const auto node3 = 2u * i;
        const auto node4 = 2u * i - 1u;
        createElement(domain, node1, node2, node3, node4, material.name());
    }

    /// Create Boundary Conditions
    const auto hc = 5.0;
    const auto airTemperature = 10.0;
    const auto humidity = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, humidity};

    createBC_FixedHc(domain, 1, 2, bcCoeff);
    createBC_FixedHc(domain, 5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 24;

    auto temperatures = HygroThermFEM::properties(HygroThermFEM::Variable::temperature);
    auto humidities = HygroThermFEM::properties(HygroThermFEM::Variable::humidity);
    std::vector<std::vector<double>> temperatureSolution;
    std::vector<std::vector<double>> waterContentSolution;

    for(auto i = 0; i < nSteps; ++i)
    {
        HygroThermFEM::TransientSubstitutionSolver solver;
        auto aSolution = solver.transient(domain, temperatures, humidities, dTime);
        temperatureSolution.push_back(aSolution.temperature);
        waterContentSolution.push_back(aSolution.waterContent);
        temperatures = aSolution.temperature;
        humidities = aSolution.humidity;
    }

    std::vector<std::vector<double>> correctWaterContentSolution{
      {13.313657, 13.313657, 25.249587, 25.249587, 13.313657, 13.313657},
      {12.614013, 12.614013, 25.243483, 25.243483, 12.614013, 12.614013},
      {11.944151, 11.944151, 25.236387, 25.236387, 11.944151, 11.944151},
      {11.503568, 11.503568, 25.228286, 25.228286, 11.503568, 11.503568},
      {11.062090, 11.062090, 25.219167, 25.219167, 11.062090, 11.062090},
      {10.620303, 10.620303, 25.209015, 25.209015, 10.620303, 10.620303},
      {10.178924, 10.178924, 25.197816, 25.197816, 10.178924, 10.178924},
      {9.742722, 9.742722, 25.185571, 25.185571, 9.742722, 9.742722},
      {9.311973, 9.311973, 25.172275, 25.172275, 9.311973, 9.311973},
      {8.886901, 8.886901, 25.157925, 25.157925, 8.886901, 8.886901},
      {8.467687, 8.467687, 25.142517, 25.142517, 8.467687, 8.467687},
      {8.102494, 8.102494, 25.126047, 25.126047, 8.102494, 8.102494},
      {7.751972, 7.751972, 25.108512, 25.108512, 7.751972, 7.751972},
      {7.406807, 7.406807, 25.089909, 25.089909, 7.406807, 7.406807},
      {7.067072, 7.067072, 25.070233, 25.070233, 7.067072, 7.067072},
      {6.732706, 6.732706, 25.049481, 25.049481, 6.732706, 6.732706},
      {6.403556, 6.403556, 25.027649, 25.027649, 6.403556, 6.403556},
      {6.079685, 6.079685, 25.004731, 25.004731, 6.079685, 6.079685},
      {5.761148, 5.761148, 24.991884, 24.991884, 5.761148, 5.761148},
      {5.447994, 5.447994, 24.981316, 24.981316, 5.447994, 5.447994},
      {5.218073, 5.218073, 24.970287, 24.970287, 5.218073, 5.218073},
      {5.063043, 5.063043, 24.958795, 24.958795, 5.063043, 5.063043},
      {4.910836, 4.910836, 24.946839, 24.946839, 4.910836, 4.910836},
      {4.761465, 4.761465, 24.934417, 24.934417, 4.761465, 4.761465},
    };

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.346950, 0.346950, 0.269697, 0.269697, 0.346950, 0.346950},
      {0.607021, 0.607021, 0.542988, 0.542988, 0.607021, 0.607021},
      {0.835412, 0.835412, 0.786306, 0.786306, 0.835412, 0.835412},
      {1.041938, 1.041938, 1.003031, 1.003031, 1.041938, 1.041938},
      {1.231312, 1.231312, 1.199297, 1.199297, 1.231312, 1.231312},
      {1.407293, 1.407293, 1.380097, 1.380097, 1.407293, 1.407293},
      {1.572976, 1.572976, 1.549261, 1.549261, 1.572976, 1.572976},
      {1.730570, 1.730570, 1.709470, 1.709470, 1.730570, 1.730570},
      {1.881848, 1.881848, 1.862763, 1.862763, 1.881848, 1.881848},
      {2.028198, 2.028198, 2.010703, 2.010703, 2.028198, 2.028198},
      {2.170701, 2.170701, 2.154487, 2.154487, 2.170701, 2.170701},
      {2.310374, 2.310374, 2.294971, 2.294971, 2.310374, 2.310374},
      {2.447516, 2.447516, 2.432998, 2.432998, 2.447516, 2.447516},
      {2.582754, 2.582754, 2.568995, 2.568995, 2.582754, 2.582754},
      {2.716471, 2.716471, 2.703367, 2.703367, 2.716471, 2.716471},
      {2.848961, 2.848961, 2.836432, 2.836432, 2.848961, 2.848961},
      {2.980461, 2.980461, 2.968441, 2.968441, 2.980461, 2.980461},
      {3.111144, 3.111144, 3.099582, 3.099582, 3.111144, 3.111144},
      {3.241133, 3.241133, 3.229985, 3.229985, 3.241133, 3.241133},
      {3.370510, 3.370510, 3.359742, 3.359742, 3.370510, 3.370510},
      {3.499313, 3.499313, 3.488896, 3.488896, 3.499313, 3.499313},
      {3.627559, 3.627559, 3.617470, 3.617470, 3.627559, 3.627559},
      {3.755271, 3.755271, 3.745488, 3.745488, 3.755271, 3.755271},
      {3.882452, 3.882452, 3.872958, 3.872958, 3.882452, 3.882452},
    };

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
