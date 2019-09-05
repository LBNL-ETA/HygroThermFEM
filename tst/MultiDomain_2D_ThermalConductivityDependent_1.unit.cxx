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
      {13.313626, 13.313626, 25.249588, 25.249588, 13.313626, 13.313626},
      {12.614009, 12.614009, 25.243483, 25.243483, 12.614009, 12.614009},
      {11.944159, 11.944159, 25.236386, 25.236386, 11.944159, 11.944159},
      {11.503585, 11.503585, 25.228286, 25.228286, 11.503585, 11.503585},
      {11.062113, 11.062113, 25.219166, 25.219166, 11.062113, 11.062113},
      {10.620331, 10.620331, 25.209013, 25.209013, 10.620331, 10.620331},
      {10.178957, 10.178957, 25.197814, 25.197814, 10.178957, 10.178957},
      {9.742758, 9.742758, 25.185568, 25.185568, 9.742758, 9.742758},
      {9.312012, 9.312012, 25.172272, 25.172272, 9.312012, 9.312012},
      {8.886943, 8.886943, 25.157922, 25.157922, 8.886943, 8.886943},
      {8.467731, 8.467731, 25.142514, 25.142514, 8.467731, 8.467731},
      {8.102506, 8.102506, 25.126044, 25.126044, 8.102506, 8.102506},
      {7.751962, 7.751962, 25.10851, 25.10851, 7.751962, 7.751962},
      {7.406778, 7.406778, 25.089907, 25.089907, 7.406778, 7.406778},
      {7.067026, 7.067026, 25.070231, 25.070231, 7.067026, 7.067026},
      {6.732646, 6.732646, 25.04948, 25.04948, 6.732646, 6.732646},
      {6.403483, 6.403483, 25.027647, 25.027647, 6.403483, 6.403483},
      {6.0796, 6.0796, 25.00473, 25.00473, 6.0796, 6.0796},
      {5.761053, 5.761053, 24.991883, 24.991883, 5.761053, 5.761053},
      {5.44789, 5.44789, 24.981316, 24.981316, 5.44789, 5.44789},
      {5.218015, 5.218015, 24.970286, 24.970286, 5.218015, 5.218015},
      {5.062982, 5.062982, 24.958794, 24.958794, 5.062982, 5.062982},
      {4.910771, 4.910771, 24.946838, 24.946838, 4.910771, 4.910771},
      {4.761397, 4.761397, 24.934417, 24.934417, 4.761397, 4.761397}};

    EXPECT_EQ(waterContentSolution.size(), correctWaterContentSolution.size());

    for(auto i = 0u; i < waterContentSolution.size(); ++i)
    {
        for(auto j = 0u; j < waterContentSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctWaterContentSolution[i][j], waterContentSolution[i][j], 1e-6);
        }
    }

    std::vector<std::vector<double>> correctTemperatureSolution{
      {0.346935, 0.346935, 0.269693, 0.269693, 0.346935, 0.346935},
      {0.607101, 0.607101, 0.542993, 0.542993, 0.607101, 0.607101},
      {0.835514, 0.835514, 0.786335, 0.786335, 0.835514, 0.835514},
      {1.042048, 1.042048, 1.003072, 1.003072, 1.042048, 1.042048},
      {1.231422, 1.231422, 1.199342, 1.199342, 1.231422, 1.231422},
      {1.4074, 1.4074, 1.380141, 1.380141, 1.4074, 1.4074},
      {1.573078, 1.573078, 1.549302, 1.549302, 1.573078, 1.573078},
      {1.730666, 1.730666, 1.709506, 1.709506, 1.730666, 1.730666},
      {1.881937, 1.881937, 1.862794, 1.862794, 1.881937, 1.881937},
      {2.02828, 2.02828, 2.010728, 2.010728, 2.02828, 2.02828},
      {2.170778, 2.170778, 2.154506, 2.154506, 2.170778, 2.170778},
      {2.310252, 2.310252, 2.295032, 2.295032, 2.310252, 2.310252},
      {2.447346, 2.447346, 2.433007, 2.433007, 2.447346, 2.447346},
      {2.582556, 2.582556, 2.568966, 2.568966, 2.582556, 2.582556},
      {2.716255, 2.716255, 2.703314, 2.703314, 2.716255, 2.716255},
      {2.848736, 2.848736, 2.836364, 2.836364, 2.848736, 2.848736},
      {2.980235, 2.980235, 2.968365, 2.968365, 2.980235, 2.980235},
      {3.110922, 3.110922, 3.099504, 3.099504, 3.110922, 3.110922},
      {3.240917, 3.240917, 3.229908, 3.229908, 3.240917, 3.240917},
      {3.370304, 3.370304, 3.359669, 3.359669, 3.370304, 3.370304},
      {3.499118, 3.499118, 3.488829, 3.488829, 3.499118, 3.499118},
      {3.627377, 3.627377, 3.617411, 3.617411, 3.627377, 3.627377},
      {3.755101, 3.755101, 3.745437, 3.745437, 3.755101, 3.755101},
      {3.882296, 3.882296, 3.872915, 3.872915, 3.882296, 3.882296}};

    EXPECT_EQ(temperatureSolution.size(), correctTemperatureSolution.size());

    for(auto i = 0u; i < temperatureSolution.size(); ++i)
    {
        for(auto j = 0u; j < temperatureSolution[i].size(); ++j)
        {
            EXPECT_NEAR(correctTemperatureSolution[i][j], temperatureSolution[i][j], 1e-6);
        }
    }
}
