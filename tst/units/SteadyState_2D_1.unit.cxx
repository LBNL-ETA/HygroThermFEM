#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

class SteadyState_2D_1 : public testing::Test
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

TEST_F(SteadyState_2D_1, TestExample_1)
{
    const double initialTemperature = 21.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;
    constexpr auto liquidPercent = 1.0;

    auto state = State(initialTemperature, initialMoistureContent, initialPressure, liquidPercent);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    NodePool::Instance().createNode(
      1, 1, 5, State(initialTemperature, 0, initialPressure, liquidPercent));
    NodePool::Instance().createNode(
      2, 1, 0, State(initialTemperature, 0, initialPressure, liquidPercent));
    NodePool::Instance().createNode(
      3, 0.5, 5, State(initialTemperature, 0.5, initialPressure, liquidPercent));
    NodePool::Instance().createNode(
      4, 0.5, 0, State(initialTemperature, 0.5, initialPressure, liquidPercent));
    NodePool::Instance().createNode(
      5, 0, 5, State(initialTemperature, 1, initialPressure, liquidPercent));
    NodePool::Instance().createNode(
      6, 0, 0, State(initialTemperature, 1, initialPressure, liquidPercent));

    // Material Properties
    constexpr double thermalConductivityDry{1.0};
    constexpr double density{2050.0};
    constexpr double porosity{0.18};
    constexpr double specificHeatCapacityDry{850.0};
    constexpr double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.0}, {180, 1.0}};
    constexpr double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.0}, {1, 1.0}};
    constexpr double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0}, {180, 7E-7}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0}, {1, 180}};

    auto & material =
      MaterialPool::Instance().createSolidMaterial("Test Material",
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

    const auto simulateThermal{true};
    const auto simulateMoisture{false};

    HygroThermFEM::MultiDomain domain{simulateThermal, simulateMoisture};

    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    // constexpr auto hc1 = 20.0;
    const auto hc1 = 1e20;
    constexpr auto humidity1 = 0.0;
    // constexpr auto temperatureAir1 = -18.0;
    constexpr auto temperatureAir1 = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1, humidity1};

    // constexpr auto hc2 = 2.4;
    const auto hc2 = 1e20;
    constexpr auto humidity2 = 1.0;
    constexpr auto temperatureAir2 = 20.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2, humidity2};

    domain.createBC_FixedHc(1, 2, bcCoeff1);
    domain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = domain.steadyState();
    const auto temperature = solution.temperature;
    const auto humidity = solution.humidity;

    std::vector<double> correctTemperature{0, 0, 10.658980, 10.658980, 20, 20};

    EXPECT_EQ(temperature.size(), correctTemperature.size());

    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(temperature[i], correctTemperature[i], 1e-6);
    }

    std::vector<double> correctHumidity{0, 0, 0.5, 0.5, 1.0, 1.0};

    EXPECT_EQ(humidity.size(), correctHumidity.size());

    for(auto i = 0u; i < correctHumidity.size(); ++i)
    {
        EXPECT_NEAR(humidity[i], correctHumidity[i], 1e-6);
    }
}
