#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::State;

TEST(SteadyState_2D_1, TestExample_1)
{
    const double initialTemperature = 21.0;
    const double initialPressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

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

    // Material Properties (using C++20 designated initializers)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Test Material",
        .thermalConductivityDry = 1.0,
        .density = 2050.0,
        .porosity = 0.18,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.0}, {180, 1.0}},
        .moistureDependentMeasurementTemperature = 0.0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.0}, {1, 1.0}},
        .temperatureDependentMeasurementHumidity = 0.0,
        .liquidTransportCurve = {{0, 0}, {180, 7E-7}},
        .sorptionCurve = {{0, 0}, {1, 180}}
    });

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

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

    multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = multiDomain.steadyState();
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
