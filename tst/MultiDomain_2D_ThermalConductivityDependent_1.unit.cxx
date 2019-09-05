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
      {0.0, 1.8}, {180, 3.0}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.2}};
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

    domain.createMoistureBCFixedHc(1, 2, bcCoeff);
    domain.createMoistureBCFixedHc(5, 6, bcCoeff);

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
      {13.310557, 13.310557, 25.250274, 25.250274, 13.310557, 13.310557},
      {12.610226, 12.610226, 25.244226, 25.244226, 12.610226, 12.610226},
      {11.941371, 11.941371, 25.237199, 25.237199, 11.941371, 11.941371},
      {11.500443, 11.500443, 25.229173, 25.229173, 11.500443, 11.500443},
      {11.058698, 11.058698, 25.220130, 25.220130, 11.058698, 11.058698},
      {10.616714, 10.616714, 25.210057, 25.210057, 10.616714, 10.616714},
      {10.175234, 10.175234, 25.198938, 25.198938, 10.175234, 10.175234},
      {9.7389750, 9.7389750, 25.186776, 25.186776, 9.7389750, 9.7389750},
      {9.308205, 9.308205, 25.173565, 25.173565, 9.308205, 9.308205},
      {8.883136, 8.883136, 25.159303, 25.159303, 8.883136, 8.883136},
      {8.463944, 8.463944, 25.143986, 25.143986, 8.463944, 8.463944},
      {8.099274, 8.099274, 25.12761, 25.12761, 8.099274, 8.099274},
      {7.748769, 7.748769, 25.110172, 25.110172, 7.748769, 7.748769},
      {7.40363, 7.40363, 25.091668, 25.091668, 7.40363, 7.40363},
      {7.063927, 7.063927, 25.072096, 25.072096, 7.063927, 7.063927},
      {6.729597, 6.729597, 25.051451, 25.051451, 6.729597, 6.729597},
      {6.400488, 6.400488, 25.029729, 25.029729, 6.400488, 6.400488},
      {6.076661, 6.076661, 25.006925, 25.006925, 6.076661, 6.076661},
      {5.758172, 5.758172, 24.992857, 24.992857, 5.758172, 5.758172},
      {5.44507, 5.44507, 24.98234, 24.98234, 5.44507, 5.44507},
      {5.2166, 5.2166, 24.971363, 24.971363, 5.2166, 5.2166},
      {5.0616, 5.0616, 24.959925, 24.959925, 5.0616, 5.0616},
      {4.909422, 4.909422, 24.948024, 24.948024, 4.909422, 4.909422},
      {4.760083, 4.760083, 24.93566, 24.93566, 4.760083, 4.760083}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.357576, 0.357576, 0.259395, 0.259395, 0.357576, 0.357576},
      {0.622865, 0.622865, 0.523134, 0.523134, 0.622865, 0.622865},
      {0.850831, 0.850831, 0.760992, 0.760992, 0.850831, 0.850831},
      {1.055077, 1.055077, 0.974521, 0.974521, 1.055077, 1.055077},
      {1.241964, 1.241964, 1.168774, 1.168774, 1.241964, 1.241964},
      {1.415773, 1.415773, 1.348246, 1.348246, 1.415773, 1.415773},
      {1.579717, 1.579717, 1.516501, 1.516501, 1.579717, 1.579717},
      {1.735961, 1.735961, 1.676091, 1.676091, 1.735961, 1.735961},
      {1.886215, 1.886215, 1.828956, 1.828956, 1.886215, 1.886215},
      {2.031805, 2.031805, 1.976591, 1.976591, 2.031805, 2.031805},
      {2.173752, 2.173752, 2.120147, 2.120147, 2.173752, 2.173752},
      {2.312828, 2.312828, 2.260500, 2.260500, 2.312828, 2.312828},
      {2.449637, 2.449637, 2.398328, 2.398328, 2.449637, 2.449637},
      {2.584643, 2.584643, 2.534156, 2.534156, 2.584643, 2.584643},
      {2.718196, 2.718196, 2.668381, 2.668381, 2.718196, 2.718196},
      {2.850572, 2.850572, 2.801313, 2.801313, 2.850572, 2.850572},
      {2.981995, 2.981995, 2.933198, 2.933198, 2.981995, 2.981995},
      {3.112624, 3.112624, 3.064220, 3.064220, 3.112624, 3.112624},
      {3.242574, 3.242574, 3.194511, 3.194511, 3.242574, 3.242574},
      {3.371922, 3.371922, 3.324164, 3.324164, 3.371922, 3.371922},
      {3.500697, 3.500697, 3.453225, 3.453225, 3.500697, 3.500697},
      {3.628916, 3.628916, 3.581718, 3.581718, 3.628916, 3.628916},
      {3.756599, 3.756599, 3.709667, 3.709667, 3.756599, 3.756599},
      {3.883751, 3.883751, 3.837081, 3.837081, 3.883751, 3.883751}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
