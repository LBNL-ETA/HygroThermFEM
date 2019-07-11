#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;
using HygroThermFEM::SimulationProperties;

class SteadyState_2D_ExcludeVaporDiffusion_1 : public testing::Test
{
protected:
    void SetUp() override
    {}

    void TearDown() override
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
        SimulationProperties::Instance().reset();
    }
};

TEST_F(SteadyState_2D_ExcludeVaporDiffusion_1, TestExample_1)
{
    const auto excludeWaterLiquidTransportation{false};
    const auto excludeHeatOfEvaporation{false};
    const auto excludeCapillaryConduction{false};
    const auto excludeVaporDiffusionConduction{true};

    SimulationProperties::Instance().setCalculationParameters(excludeWaterLiquidTransportation,
                                                              excludeHeatOfEvaporation,
                                                              excludeCapillaryConduction,
                                                              excludeVaporDiffusionConduction);

    const double initialTemperature = 21.0;
    const double initialMoistureContent = 0.0;
    const double initialPressure = 101325;
    const auto liquidPercent = 1.0;

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

    auto & material = MaterialPool::Instance().createSolidMaterial(
      "Test Material",
      2050,                       /// Density
      0.18,                       /// Porosity
      850,                        /// Specific Heat Capacity (dry)
      15E-6,                      /// Diffusion Resistance Factor
      {{0.0, 1.0}, {180, 1.0}},   /// Thermal Conductivity (as function of water content)
      {{0, 0},   /// Liquid Transportation Coefficient (as function of water content)
       {27, 1E-8},
       {45, 1.1E-8},
       {90, 2E-8},
       {126, 3.5E-8},
       {144, 5E-8},
       {162, 1E-7},
       {171, 2E-7},
       {180, 7E-7}},
      {{0, 0},   /// Moisture Storage Function (Water content as function of relative humidity)
       {0.5, 5.3},
       {0.65, 8.4},
       {0.8, 12},
       {0.93, 17},
       {0.95, 25},
       {0.99, 63},
       {0.995, 83},
       {0.999, 120},
       {1, 180}});

    const auto simulateThermal{true};
    const auto simulateMoisture{true};

    HygroThermFEM::MultiDomain domain{simulateThermal, simulateMoisture};

    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto hc1 = 1e20;
    const auto humidity1 = 0.8;
    const auto temperatureAir1 = 0.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{temperatureAir1, hc1, humidity1};

    const auto hc2 = 1e20;
    const auto humidity2 = 0.0;
    const auto temperatureAir2 = 20.0;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{temperatureAir2, hc2, humidity2};

    domain.createMoistureBCFixedHc(1, 2, bcCoeff1);
    domain.createMoistureBCFixedHc(6, 5, bcCoeff2);

    const auto solution = domain.steadyState();
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
