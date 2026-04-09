#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "ObserveSimulationProgress.hxx"
#include "TestHelpers.hxx"
#include "TestMaterials.hxx"

class MultiDomain_BeamWall_SingleMaterial : public testing::Test
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

TEST_F(MultiDomain_BeamWall_SingleMaterial, Stucco_99dot9_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.999;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the stucco material so it is available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of stucco, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, exteriorBc);
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, interiorBc);

    const auto results = multiDomain.transientMultiStep(3600.0, 10);

    const std::vector<std::vector<double>> correctWaterContentSolution{
      {110.001861, 110.001861, 110.000106, 110.000106, 110.000106, 110.000106, 110.001861, 110.001861},
      {110.002631, 110.002631, 110.000760, 110.000760, 110.000760, 110.000760, 110.002631, 110.002631},
      {110.091966, 110.091966, 110.089115, 110.089115, 110.089115, 110.089115, 110.091966, 110.091966},
      {110.092130, 110.092130, 110.091111, 110.091111, 110.091111, 110.091111, 110.092130, 110.092130},
      {110.093352, 110.093352, 110.091536, 110.091536, 110.091536, 110.091536, 110.093352, 110.093352},
      {110.094078, 110.094078, 110.092210, 110.092210, 110.092210, 110.092210, 110.094078, 110.094078},
      {110.182876, 110.182876, 110.180587, 110.180587, 110.180587, 110.180587, 110.182876, 110.182876},
      {110.184579, 110.184579, 110.183880, 110.183880, 110.183880, 110.183880, 110.184579, 110.184579},
      {110.185996, 110.185996, 110.184205, 110.184205, 110.184205, 110.184205, 110.185996, 110.185996},
      {110.186732, 110.186732, 110.184873, 110.184873, 110.184873, 110.184873, 110.186732, 110.186732}};

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{1.532058e-12, 1.485851e-12, 1.093719e-12, 1.816684e-12, 1.501880e-12, 1.481448e-12, 1.309664e-12, 1.938170e-12, 1.506435e-12, 1.479089e-12};
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}

TEST_F(MultiDomain_BeamWall_SingleMaterial, Stucco_90_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.9;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the stucco material so it is available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of stucco, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

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
      {77.708163, 77.708163, 67.058496, 67.058496, 67.058496, 67.058496, 77.708163, 77.708163},
      {81.554138, 81.554138, 71.108940, 71.108940, 71.108940, 71.108940, 81.554138, 81.554138},
      {82.808884, 82.808884, 75.875683, 75.875683, 75.875683, 75.875683, 82.808884, 82.808884},
      {83.574311, 83.574311, 80.226113, 80.226113, 80.226113, 80.226113, 83.574311, 83.574311},
      {84.853557, 84.853557, 82.496868, 82.496868, 82.496868, 82.496868, 84.853557, 84.853557},
      {86.019139, 86.019139, 84.103202, 84.103202, 84.103202, 84.103202, 86.019139, 86.019139},
      {87.017049, 87.017049, 85.371609, 85.371609, 85.371609, 85.371609, 87.017049, 87.017049},
      {87.876194, 87.876194, 86.426443, 86.426443, 86.426443, 86.426443, 87.876194, 87.876194},
      {88.624896, 88.624896, 87.330480, 87.330480, 87.330480, 87.330480, 88.624896, 88.624896},
      {89.286738, 89.286738, 88.119948, 88.119948, 88.119948, 88.119948, 89.286738, 89.286738}};

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{2.431356e-08, 1.267140e-08, 9.515974e-09, 8.041952e-09, 5.562765e-09, 4.323681e-09, 3.527430e-09, 2.969864e-09, 2.558945e-09, 2.241831e-09};
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}

TEST_F(MultiDomain_BeamWall_SingleMaterial, Stucco_50_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.5;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the stucco material so it is available to the elements.
    const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of stucco, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = stucco.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

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
      {46.082213, 46.082213, 31.436613, 31.436613, 31.436613, 31.436613, 46.082213, 46.082213},
      {63.248663, 63.248663, 33.888978, 33.888978, 33.888978, 33.888978, 63.248663, 63.248663},
      {72.749372, 72.749372, 37.036583, 37.036583, 37.036583, 37.036583, 72.749372, 72.749372},
      {77.346099, 77.346099, 39.516176, 39.516176, 39.516176, 39.516176, 77.346099, 77.346099},
      {79.648597, 79.648597, 41.510678, 41.510678, 41.510678, 41.510678, 79.648597, 79.648597},
      {80.446839, 80.446839, 43.276153, 43.276153, 43.276153, 43.276153, 80.446839, 80.446839},
      {80.795132, 80.795132, 44.856190, 44.856190, 44.856190, 44.856190, 80.795132, 80.795132},
      {81.090558, 81.090558, 50.103310, 50.103310, 50.103310, 50.103310, 81.090558, 81.090558},
      {81.390458, 81.390458, 55.254292, 55.254292, 55.254292, 55.254292, 81.390458, 81.390458},
      {81.734477, 81.734477, 59.727555, 59.727555, 59.727555, 59.727555, 81.734477, 81.734477}};

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{2.813586e-07, 9.292696e-08, 5.533571e-08, 3.481325e-08, 2.502997e-08, 1.977597e-08, 1.707836e-08, 1.536177e-08, 1.401509e-08, 1.231617e-08};
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}

TEST_F(MultiDomain_BeamWall_SingleMaterial, Fiberglass_99dot99999_Percent)
{
    constexpr auto excludeWaterLiquidTransportation{false};
    constexpr auto excludeHeatOfEvaporation{true};
    constexpr auto excludeCapillaryConduction{true};
    constexpr auto excludeVaporDiffusionConduction{true};
    constexpr auto thermalConductivityMoistureAndTemperatureDependent{false};
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      excludeWaterLiquidTransportation,
      excludeHeatOfEvaporation,
      excludeCapillaryConduction,
      excludeVaporDiffusionConduction,
      thermalConductivityMoistureAndTemperatureDependent);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 30.0;
    constexpr double initialDomainHumidity = 0.9999999;
    constexpr double bcHumidity = 1.0;

    // Initial state: warm, near-saturation (the regime where the solver
    // currently struggles).
    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    // Register the fiberglass material so it is available to the elements.
    const auto & fiberglass =
      multiDomain.materials().createSolidMaterial(TestHelper::FiberglassBatts());

    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);

    // Beam geometry:
    //   - 0.02 m of fiberglass, 3 elements wide
    //   - Total height 0.05 m, 1 element row
    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = fiberglass.name(), .numElementsX = 3, .width = 0.02})
      .build();

    // Boundary conditions: identical to the initial domain state on both
    // edges. With no thermal or moisture gradient driving the problem the
    // solution should remain constant; any drift exposes solver issues.
    constexpr auto hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{
      /*airTemperature=*/initialTemperature, hc, /*airHumidity=*/bcHumidity};

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
      {13.399900, 13.399900, 13.399881, 13.399881, 13.399881, 13.399881, 13.399900, 13.399900},
      {13.399915, 13.399915, 13.399882, 13.399882, 13.399882, 13.399882, 13.399915, 13.399915},
      {13.399927, 13.399927, 13.399884, 13.399884, 13.399884, 13.399884, 13.399927, 13.399927},
      {13.399936, 13.399936, 13.399886, 13.399886, 13.399886, 13.399886, 13.399936, 13.399936},
      {13.399943, 13.399943, 13.399888, 13.399888, 13.399888, 13.399888, 13.399943, 13.399943},
      {13.399949, 13.399949, 13.399891, 13.399891, 13.399891, 13.399891, 13.399949, 13.399949},
      {13.399954, 13.399954, 13.399893, 13.399893, 13.399893, 13.399893, 13.399954, 13.399954},
      {13.399957, 13.399957, 13.399896, 13.399896, 13.399896, 13.399896, 13.399957, 13.399957},
      {13.399961, 13.399961, 13.399899, 13.399899, 13.399899, 13.399899, 13.399961, 13.399961},
      {13.399963, 13.399963, 13.399901, 13.399901, 13.399901, 13.399901, 13.399963, 13.399963}};

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector correctHumidityError(10u, 0.0);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}

TEST_F(MultiDomain_BeamWall_SingleMaterial, GypsumBoardHumidityGradient)
{
    HygroThermFEM::SimulationProperties::Instance().setCalculationParameters(
      /*excludeWaterLiquidTransportation=*/false,
      /*excludeHeatOfEvaporation=*/false,
      /*excludeCapillaryConduction=*/false,
      /*excludeVaporDiffusionConduction=*/false,
      /*thermalConductivityMoistureAndTemperatureDependent=*/false);

    HygroThermFEM::MultiDomain multiDomain({.performThermal = false, .performMoisture = true});

    constexpr double initialTemperature = 20.0;
    constexpr double initialDomainHumidity = 0.99999;
    constexpr double bcHumidity = 1.0;

    constexpr HygroThermFEM::State initialState({.temperature = initialTemperature,
                                                 .humidity = initialDomainHumidity,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

    const auto & material = multiDomain.materials().createSolidMaterial(TestHelper::GypsumBoardInterior());

    TestHelper::ObserveSimulationProgress progressMoisture;
    multiDomain.subscribeMoisture(&progressMoisture);
    TestHelper::ObserveSimulationProgress progressThermal;
    multiDomain.subscribeThermal(&progressThermal);

    TestHelper::BeamBuilder builder(multiDomain);
    builder.xStart(0.0)
      .height(0.05)
      .numElementsY(1)
      .state(initialState)
      .addSegment({.material = material.name(), .numElementsX = 5, .width = 0.10})
      .build();

    constexpr double hc = 10.0;
    const HygroThermFEM::FixedBCHCCoefficients exteriorBc{/*airTemperature=*/initialTemperature,
                                                          hc,
                                                          /*airHumidity=*/bcHumidity};
    const HygroThermFEM::FixedBCHCCoefficients interiorBc{/*airTemperature=*/initialTemperature,
                                                          hc,
                                                          /*airHumidity=*/bcHumidity};

    for(auto [i1, i2] : builder.leftEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, exteriorBc);
    }
    for(auto [i1, i2] : builder.rightEdge())
    {
        multiDomain.createBC_FixedHc(i1, i2, interiorBc);
    }

    constexpr double dTime = 3600.0;
    constexpr int nSteps = 20;

    const auto results = multiDomain.transientMultiStep(dTime, nSteps);

    const std::vector<std::vector<double>> correctWaterContentSolution{
      {355.000002, 355.000002, 355.000001, 355.000001, 355.000000, 355.000000, 355.000000, 355.000000, 355.000001, 355.000001, 355.000002, 355.000002},
      {354.999987, 354.999987, 355.000019, 355.000019, 355.000045, 355.000045, 355.000045, 355.000045, 355.000019, 355.000019, 354.999987, 354.999987},
      {355.008680, 355.008680, 354.999864, 354.999864, 354.999654, 354.999654, 354.999654, 354.999654, 354.999864, 354.999864, 355.008680, 355.008680},
      {366.117051, 366.117051, 366.163208, 366.163208, 366.183604, 366.183604, 366.183604, 366.183604, 366.163208, 366.163208, 366.117051, 366.117051},
      {369.969917, 369.969917, 369.970005, 369.970005, 369.970079, 369.970079, 369.970079, 369.970079, 369.970005, 369.970005, 369.969917, 369.969917},
      {369.970152, 369.970152, 369.970239, 369.970239, 369.965666, 369.965666, 369.965666, 369.965666, 369.970239, 369.970239, 369.970152, 369.970152},
      {369.864120, 369.864120, 369.817558, 369.817558, 369.965934, 369.965934, 369.965934, 369.965934, 369.817558, 369.817558, 369.864120, 369.864120},
      {369.998978, 369.998978, 369.998992, 369.998992, 369.999297, 369.999297, 369.999297, 369.999297, 369.998992, 369.998992, 369.998978, 369.998978},
      {369.998986, 369.998986, 369.999000, 369.999000, 369.985323, 369.985323, 369.985323, 369.985323, 369.999000, 369.999000, 369.998986, 369.998986},
      {337.534479, 337.534479, 333.561719, 333.561719, 360.946547, 360.946547, 360.946547, 360.946547, 333.561719, 333.561719, 337.534479, 337.534479},
      {337.089301, 337.089301, 332.338215, 332.338215, 357.909227, 357.909227, 357.909227, 357.909227, 332.338215, 332.338215, 337.089301, 337.089301},
      {336.955793, 336.955793, 332.185415, 332.185415, 357.308964, 357.308964, 357.308964, 357.308964, 332.185415, 332.185415, 336.955793, 336.955793},
      {336.941692, 336.941692, 332.163499, 332.163499, 357.216818, 357.216818, 357.216818, 357.216818, 332.163499, 332.163499, 336.941692, 336.941692},
      {336.939834, 336.939834, 332.161028, 332.161028, 357.205716, 357.205716, 357.205716, 357.205716, 332.161028, 332.161028, 336.939834, 336.939834},
      {336.939628, 336.939628, 332.160733, 332.160733, 357.204403, 357.204403, 357.204403, 357.204403, 332.160733, 332.160733, 336.939628, 336.939628},
      {336.939604, 336.939604, 332.160699, 332.160699, 357.204251, 357.204251, 357.204251, 357.204251, 332.160699, 332.160699, 336.939604, 336.939604},
      {336.939601, 336.939601, 332.160695, 332.160695, 357.204233, 357.204233, 357.204233, 357.204233, 332.160695, 332.160695, 336.939601, 336.939601},
      {336.939601, 336.939601, 332.160694, 332.160694, 357.204231, 357.204231, 357.204231, 357.204231, 332.160694, 332.160694, 336.939601, 336.939601},
      {336.939600, 336.939600, 332.160694, 332.160694, 357.204231, 357.204231, 357.204231, 357.204231, 332.160694, 332.160694, 336.939600, 336.939600},
      {336.939600, 336.939600, 332.160694, 332.160694, 357.204231, 357.204231, 357.204231, 357.204231, 332.160694, 332.160694, 336.939600, 336.939600}};

    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{0.000000, 4.338247e-13, 5.294456e-11, 1.919323e-10, 3.927972e-13, 3.897204e-13, 2.075163e-10, 6.179120e-14, 6.115021e-14, 8.919306e-07, 8.265528e-07, 8.138665e-07, 8.118773e-07, 8.116380e-07, 8.116096e-07, 8.116063e-07, 8.116059e-07, 8.116059e-07, 8.116059e-07, 8.116059e-07};
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(20u, std::vector(12u, 20.0));
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(20u, 0.0);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
