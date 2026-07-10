#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::State;
using HygroThermFEM::SimulationProperties;

class SteadyState_2D_ExcludeLiquidTransport_1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(SteadyState_2D_ExcludeLiquidTransport_1, TestExample_1)
{
    constexpr auto excludeWaterLiquidTransportation{true};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    constexpr double initialTemperature = 21.0;
    constexpr double initialPressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    multiDomain.nodes().createNode({.index = 1, .x = 1, .y = 5, .state = State{
        .temperature = initialTemperature,
        .humidity = 0,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }});
    multiDomain.nodes().createNode({.index = 2, .x = 1, .y = 0, .state = State{
        .temperature = initialTemperature,
        .humidity = 0,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }});
    multiDomain.nodes().createNode({.index = 3, .x = 0.5, .y = 5, .state = State{
        .temperature = initialTemperature,
        .humidity = 0.5,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }});
    multiDomain.nodes().createNode({.index = 4, .x = 0.5, .y = 0, .state = State{
        .temperature = initialTemperature,
        .humidity = 0.5,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5, .state = State{
        .temperature = initialTemperature,
        .humidity = 1,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = State{
        .temperature = initialTemperature,
        .humidity = 1,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }});

    auto params = TestHelper::TestMaterial();
    params.porosity = 0.18;
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Create Boundary Conditions
    // constexpr auto hc1 = 20.0;
    constexpr auto hc1 = 1e20;
    constexpr auto humidity1 = 0.8;
    // constexpr auto temperatureAir1 = -18.0;
    constexpr auto temperatureAir1 = 0.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1, humidity1};

    // constexpr auto hc2 = 2.4;
    constexpr auto hc2 = 1e20;
    constexpr auto humidity2 = 0.0;
    constexpr auto temperatureAir2 = 20.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2, humidity2};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = multiDomain.steadyState();
    const auto temperature = solution.temperature;
    const auto humidity = solution.humidity;

    std::vector<double> correctTemperature{-1.67840876e-15, 8.39504378e-16, 9.99999513, 9.99999513, 20, 20};

    TestHelper::dumpGolden("correctTemperature", temperature);
    TestHelper::dumpGolden("correctHumidity", humidity);
    EXPECT_EQ(temperature.size(), correctTemperature.size());

    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(temperature[i], correctTemperature[i], 1e-6);
    }

    std::vector<double> correctHumidity{0.8, 0.8, 0.201468734, 0.201468734, 4.66016607e-24, 4.66016607e-24};

    EXPECT_EQ(humidity.size(), correctHumidity.size());

    for(auto i = 0u; i < correctHumidity.size(); ++i)
    {
        EXPECT_NEAR(humidity[i], correctHumidity[i], 1e-6);
    }
}
