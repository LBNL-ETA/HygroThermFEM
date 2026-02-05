#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(ConvectionBC_2D_SteadyState, TestExample_1)
{
    SCOPED_TRACE("Begin Test: Two elementsCreator example with simple conduction.");

    HygroThermFEM::MultiDomain multiDomain;

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate

    multiDomain.nodes().createNode({.x = 15, .y = 5});
    multiDomain.nodes().createNode({.x = 15, .y = 0});
    multiDomain.nodes().createNode({.x = 5, .y = 5});
    multiDomain.nodes().createNode({.x = 5, .y = 0});
    multiDomain.nodes().createNode({.x = 0, .y = 5});
    multiDomain.nodes().createNode({.x = 0, .y = 0});

    // Material Properties (using C++20 designated initializers)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Test Material",
        .thermalConductivityDry = 1.0,
        .density = 2050.0,
        .porosity = 0.0,
        .heatCapacity = 850.0,
        .diffusionResistanceFactor = 15.0,
        .thermalConductivityMoistureDependent = {{0.0, 1.0}, {180, 1.0}},
        .moistureDependentMeasurementTemperature = 0.0,
        .thermalConductivityTemperatureDependent = {{0.0, 1.0}, {1, 1.0}},
        .temperatureDependentMeasurementHumidity = 0.0,
        .liquidTransportCurve = {{0, 0}, {27, 1E-8}, {45, 1.1E-8}, {90, 2E-8}, {126, 3.5E-8},
                                 {144, 5E-8}, {162, 1E-7}, {171, 2E-7}, {180, 7E-7}},
        .sorptionCurve = {{0, 0}, {0.5, 5.3}, {0.65, 8.4}, {0.8, 12}, {0.93, 17},
                          {0.95, 25}, {0.99, 63}, {0.995, 83}, {0.999, 120}, {1, 180}}
    });

    multiDomain.createElement(3, 4, 2, 1, material.name());
    multiDomain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    constexpr auto hc1 = 20.0;
    constexpr auto humidity1 = 0.0;
    constexpr auto temperatureAir1 = -18.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1, humidity1};

    constexpr auto hc2 = 2.4;
    constexpr auto humidity2 = 0.0;
    constexpr auto temperatureAir2 = 21.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2, humidity2};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = multiDomain.steadyState();
    auto temperature = solution.temperature;

    std::vector<double> correctTemperature{
      -17.87392241, -17.87392241, 7.341594828, 7.341594828, 19.94935345, 19.94935345};

    EXPECT_EQ(temperature.size(), correctTemperature.size());

    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(temperature[i], correctTemperature[i], 1e-6);
    }
}
