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

    const std::vector<std::vector<double>> correctWaterContentSolution{{110.227074, 110.227074, 110.012395, 110.012395, 110.012395, 110.012395, 110.227074, 110.227074},
 {110.432115, 110.432115, 110.035373, 110.035373, 110.035373, 110.035373, 110.432115, 110.432115},
 {110.618332, 110.618332, 110.067372, 110.067372, 110.067372, 110.067372, 110.618332, 110.618332},
 {110.78847, 110.78847, 110.107052, 110.107052, 110.107052, 110.107052, 110.78847, 110.78847},
 {110.944873, 110.944873, 110.15327, 110.15327, 110.15327, 110.15327, 110.944873, 110.944873},
 {111.089553, 111.089553, 110.205044, 110.205044, 110.205044, 110.205044, 111.089553, 111.089553},
 {111.224231, 111.224231, 110.261535, 110.261535, 110.261535, 110.261535, 111.224231, 111.224231},
 {111.350379, 111.350379, 110.322027, 110.322027, 110.322027, 110.322027, 111.350379, 111.350379},
 {111.469252, 111.469252, 110.385904, 110.385904, 110.385904, 110.385904, 111.469252, 111.469252},
 {111.581923, 111.581923, 110.452646, 110.452646, 110.452646, 110.452646, 111.581923, 111.581923}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{1.99742394e-12, 1.90217767e-12, 1.82044814e-12, 1.75034926e-12, 1.69015212e-12, 1.63875648e-12, 1.59474782e-12, 1.55686876e-12, 1.52464777e-12, 1.49651318e-12};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{78.8557341, 78.8557341, 65.458616, 65.458616, 65.458616, 65.458616, 78.8557341, 78.8557341},
 {85.8867038, 85.8867038, 66.6632621, 66.6632621, 66.6632621, 66.6632621, 85.8867038, 85.8867038},
 {89.2945202, 89.2945202, 68.3871397, 68.3871397, 68.3871397, 68.3871397, 89.2945202, 89.2945202},
 {90.9897313, 90.9897313, 70.3244458, 70.3244458, 70.3244458, 70.3244458, 90.9897313, 90.9897313},
 {91.9127677, 91.9127677, 72.2977964, 72.2977964, 72.2977964, 72.2977964, 91.9127677, 91.9127677},
 {92.4927591, 92.4927591, 74.2227196, 74.2227196, 74.2227196, 74.2227196, 92.4927591, 92.4927591},
 {92.9193118, 92.9193118, 76.0626007, 76.0626007, 76.0626007, 76.0626007, 92.9193118, 92.9193118},
 {93.2750176, 93.2750176, 77.8030111, 77.8030111, 77.8030111, 77.8030111, 93.2750176, 93.2750176},
 {93.5955173, 93.5955173, 79.4394816, 79.4394816, 79.4394816, 79.4394816, 93.5955173, 93.5955173},
 {93.856273, 93.856273, 81.0069377, 81.0069377, 81.0069377, 81.0069377, 93.856273, 93.856273}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{2.37656229e-08, 1.35049518e-08, 8.255845e-09, 5.73520859e-09, 4.51138601e-09, 3.87032394e-09, 3.48615653e-09, 3.2156643e-09, 2.9968247e-09, 2.79356197e-09};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{58.9855178, 58.9855178, 31.9496404, 31.9496404, 31.9496404, 31.9496404, 58.9855178, 58.9855178},
 {73.0272196, 73.0272196, 34.5256109, 34.5256109, 34.5256109, 34.5256109, 73.0272196, 73.0272196},
 {80.635428, 80.635428, 37.4330771, 37.4330771, 37.4330771, 37.4330771, 80.635428, 80.635428},
 {84.1318537, 84.1318537, 41.0704863, 41.0704863, 41.0704863, 41.0704863, 84.1318537, 84.1318537},
 {86.3661104, 86.3661104, 44.4916836, 44.4916836, 44.4916836, 44.4916836, 86.3661104, 86.3661104},
 {87.6515718, 87.6515718, 47.8997799, 47.8997799, 47.8997799, 47.8997799, 87.6515718, 87.6515718},
 {88.4645679, 88.4645679, 51.2358103, 51.2358103, 51.2358103, 51.2358103, 88.4645679, 88.4645679},
 {89.1002601, 89.1002601, 54.4194192, 54.4194192, 54.4194192, 54.4194192, 89.1002601, 89.1002601},
 {89.664783, 89.664783, 57.424529, 57.424529, 57.424529, 57.424529, 89.664783, 89.664783},
 {90.2038175, 90.2038175, 60.2379652, 60.2379652, 60.2379652, 60.2379652, 90.2038175, 90.2038175}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.27347452e-07, 6.69078963e-08, 4.50374683e-08, 4.35086447e-08, 3.89524409e-08, 1.48859468e-08, 9.94225949e-09, 9.22364622e-09, 8.6312543e-09, 8.08546871e-09};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{13.3998995, 13.3998995, 13.3998808, 13.3998808, 13.3998808, 13.3998808, 13.3998995, 13.3998995},
 {13.3999146, 13.3999146, 13.3998823, 13.3998823, 13.3998823, 13.3998823, 13.3999146, 13.3999146},
 {13.3999264, 13.3999264, 13.3998842, 13.3998842, 13.3998842, 13.3998842, 13.3999264, 13.3999264},
 {13.3999356, 13.3999356, 13.3998863, 13.3998863, 13.3998863, 13.3998863, 13.3999356, 13.3999356},
 {13.3999428, 13.3999428, 13.3998888, 13.3998888, 13.3998888, 13.3998888, 13.3999428, 13.3999428},
 {13.3999485, 13.3999485, 13.3998913, 13.3998913, 13.3998913, 13.3998913, 13.3999485, 13.3999485},
 {13.3999531, 13.3999531, 13.3998939, 13.3998939, 13.3998939, 13.3998939, 13.3999531, 13.3999531},
 {13.3999569, 13.3999569, 13.3998966, 13.3998966, 13.3998966, 13.3998966, 13.3999569, 13.3999569},
 {13.3999599, 13.3999599, 13.3998993, 13.3998993, 13.3998993, 13.3998993, 13.3999599, 13.3999599},
 {13.3999624, 13.3999624, 13.399902, 13.399902, 13.399902, 13.399902, 13.3999624, 13.3999624}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector correctHumidityError(10u, 0.0);
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(10u, std::vector(8u, 30.0));
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(10u, 0.0);
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{355.00022, 355.00022, 355.000087, 355.000087, 355.000043, 355.000043, 355.000043, 355.000043, 355.000087, 355.000087, 355.00022, 355.00022},
 {355.000359, 355.000359, 355.000186, 355.000186, 355.000114, 355.000114, 355.000114, 355.000114, 355.000186, 355.000186, 355.000359, 355.000359},
 {355.000472, 355.000472, 355.000285, 355.000285, 355.000199, 355.000199, 355.000199, 355.000199, 355.000285, 355.000285, 355.000472, 355.000472},
 {355.000575, 355.000575, 355.000383, 355.000383, 355.00029, 355.00029, 355.00029, 355.00029, 355.000383, 355.000383, 355.000575, 355.000575},
 {355.000673, 355.000673, 355.00048, 355.00048, 355.000384, 355.000384, 355.000384, 355.000384, 355.00048, 355.00048, 355.000673, 355.000673},
 {355.00077, 355.00077, 355.000576, 355.000576, 355.000479, 355.000479, 355.000479, 355.000479, 355.000576, 355.000576, 355.00077, 355.00077},
 {355.000867, 355.000867, 355.000672, 355.000672, 355.000575, 355.000575, 355.000575, 355.000575, 355.000672, 355.000672, 355.000867, 355.000867},
 {355.000963, 355.000963, 355.000768, 355.000768, 355.000671, 355.000671, 355.000671, 355.000671, 355.000768, 355.000768, 355.000963, 355.000963},
 {355.001059, 355.001059, 355.000864, 355.000864, 355.000767, 355.000767, 355.000767, 355.000767, 355.000864, 355.000864, 355.001059, 355.001059},
 {355.001155, 355.001155, 355.00096, 355.00096, 355.000863, 355.000863, 355.000863, 355.000863, 355.00096, 355.00096, 355.001155, 355.001155},
 {355.001251, 355.001251, 355.001056, 355.001056, 355.000959, 355.000959, 355.000959, 355.000959, 355.001056, 355.001056, 355.001251, 355.001251},
 {355.001347, 355.001347, 355.001152, 355.001152, 355.001055, 355.001055, 355.001055, 355.001055, 355.001152, 355.001152, 355.001347, 355.001347},
 {355.001443, 355.001443, 355.001248, 355.001248, 355.001151, 355.001151, 355.001151, 355.001151, 355.001248, 355.001248, 355.001443, 355.001443},
 {355.001539, 355.001539, 355.001344, 355.001344, 355.001247, 355.001247, 355.001247, 355.001247, 355.001344, 355.001344, 355.001539, 355.001539},
 {355.001635, 355.001635, 355.00144, 355.00144, 355.001343, 355.001343, 355.001343, 355.001343, 355.00144, 355.00144, 355.001635, 355.001635},
 {355.001732, 355.001732, 355.001536, 355.001536, 355.001439, 355.001439, 355.001439, 355.001439, 355.001536, 355.001536, 355.001732, 355.001732},
 {355.001828, 355.001828, 355.001632, 355.001632, 355.001535, 355.001535, 355.001535, 355.001535, 355.001632, 355.001632, 355.001828, 355.001828},
 {355.001924, 355.001924, 355.001728, 355.001728, 355.001631, 355.001631, 355.001631, 355.001631, 355.001728, 355.001728, 355.001924, 355.001924},
 {355.00202, 355.00202, 355.001824, 355.001824, 355.001727, 355.001727, 355.001727, 355.001727, 355.001824, 355.001824, 355.00202, 355.00202},
 {355.002116, 355.002116, 355.00192, 355.00192, 355.001823, 355.001823, 355.001823, 355.001823, 355.00192, 355.00192, 355.002116, 355.002116}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    TestHelper::dumpGolden("correctHumidityError", results.moisture.errors);
    TestHelper::expectNear(correctHumidityError, results.moisture.errors, 1e-6);

    const std::vector correctTemperatureSolution(20u, std::vector(12u, 20.0));
    TestHelper::dumpGolden("correctTemperatureSolution", results.temperature.values);
    TestHelper::expectNear(correctTemperatureSolution, results.temperature.values, 1e-6);

    const std::vector correctTemperatureError(20u, 0.0);
    TestHelper::dumpGolden("correctTemperatureError", results.temperature.errors);
    TestHelper::expectNear(correctTemperatureError, results.temperature.errors, 1e-6);

    EXPECT_EQ(progressMoisture.calls(), (std::vector<unsigned>{}));
    EXPECT_EQ(progressThermal.calls(), (std::vector<unsigned>{}));
}
