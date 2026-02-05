#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::State;
using HygroThermFEM::StateParams;
using HygroThermFEM::SimulationProperties;

class SteadyState_2D_ExcludeVaporDiffusion_1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

TEST_F(SteadyState_2D_ExcludeVaporDiffusion_1, TestExample_1)
{
    const auto excludeWaterLiquidTransportation{false};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{true};
    const auto thermalConductivityMoistureAndTemperatureDependent{false};

    SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    const double initialTemperature = 21.0;
    const double initialPressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::MultiDomain multiDomain(true, true);

    // Enter nodes. Arguments are: node number, x-coordinate, y-coordinate
    multiDomain.nodes().createNode(1, 1, 5, State({
        .temperature = initialTemperature,
        .humidity = 0,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }));
    multiDomain.nodes().createNode(2, 1, 0, State({
        .temperature = initialTemperature,
        .humidity = 0,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }));
    multiDomain.nodes().createNode(3, 0.5, 5, State({
        .temperature = initialTemperature,
        .humidity = 0.5,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }));
    multiDomain.nodes().createNode(4, 0.5, 0, State({
        .temperature = initialTemperature,
        .humidity = 0.5,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }));
    multiDomain.nodes().createNode(5, 0, 5, State({
        .temperature = initialTemperature,
        .humidity = 1,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }));
    multiDomain.nodes().createNode(6, 0, 0, State({
        .temperature = initialTemperature,
        .humidity = 1,
        .pressure = initialPressure,
        .liquidPercent = liquidPercent
    }));

    // Material Properties (Cottaer Sandstone, using C++20 designated initializers)
    const auto & material = multiDomain.materials().createSolidMaterial({
        .name = "Test material",
        .thermalConductivityDry = 1.0,
        .density = 2050.0,
        .porosity = 0.18,
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
    const auto hc1 = 1e20;
    constexpr auto humidity1 = 0.8;
    constexpr auto temperatureAir1 = 0.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1, humidity1};

    const auto hc2 = 1e20;
    constexpr auto humidity2 = 0.0;
    constexpr auto temperatureAir2 = 20.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2, humidity2};

    multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = multiDomain.steadyState();
    const auto temperature = solution.temperature;
    const auto humidity = solution.humidity;

    std::vector<double> correctTemperature{0, 0, 10, 10, 20, 20};

    EXPECT_EQ(temperature.size(), correctTemperature.size());

    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(temperature[i], correctTemperature[i], 1e-6);
    }

    std::vector<double> correctHumidity{0.8, 0.8, 0.438979, 0.438979, 0, 0};

    EXPECT_EQ(humidity.size(), correctHumidity.size());

    for(auto i = 0u; i < correctHumidity.size(); ++i)
    {
        EXPECT_NEAR(humidity[i], correctHumidity[i], 1e-6);
    }
}
