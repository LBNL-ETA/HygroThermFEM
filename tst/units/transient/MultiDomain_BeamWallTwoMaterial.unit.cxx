#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "ObserveSimulationProgress.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

class MultiDomain_BeamWall_StuccoFiberglass : public testing::Test
{
protected:
    void TearDown() override
    {
        // Reset the SimulationProperties singleton so any flags this test
        // toggled (e.g. excludeWaterLiquidTransportation) do not bleed into
        // subsequent tests in the same binary.
        HygroThermFEM::SimulationProperties::Instance().reset();
    }
};

//! Beam-shaped wall section with two materials: an exterior stucco layer
//! and a fiberglass-batt insulation layer. Initialised at near-saturation
//! humidity (RH = 0.9999) and warm temperature, with boundary conditions
//! matched exactly to the initial state. With no thermal or moisture
//! gradient driving the problem, the moisture solver must keep the
//! solution bit-stationary; any drift from the equilibrium values exposes
//! a solver- or assembly-side defect at the multi-material interface.
//! This is the regression target for the round-off filter / per-component
//! convergence metric / per-material liquid-transport coefficient fixes
//! introduced 2026-04-08.
TEST_F(MultiDomain_BeamWall_StuccoFiberglass, TwoMaterialHighHumidity)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{false};
    constexpr auto excludeCapillaryConduction{false};
    constexpr auto excludeVaporDiffusionConduction{false};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialHumidity = 0.9999;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register both materials so they are available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());
    const auto & fiberglass =
      multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of exterior stucco (left), 4 elements wide
    //   - 0.10 m of fiberglass batt insulation (right), 10 elements wide
    //   - Total height 0.05 m, 4 element rows
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .addSegment({.material = fiberglass.name(), .numElementsX = 2, .width = 0.10})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/initialHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/initialHumidity};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 10;

    const auto results = multiDomain.transientMultiStep(dTime, nSteps);

    const std::vector<std::vector<double>> correctWaterContentSolution{
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000},
      {200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 200.000000, 106.640000, 106.640000, 13.280000, 13.280000, 13.280000, 13.280000}};
    const std::vector correctHumidityError(10u, 0.0);

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector<std::vector<double>> correctTemperatureSolution{
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000},
      {30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000, 30.000000}};
    const std::vector correctTemperatureError(10u, 0.0);

    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
