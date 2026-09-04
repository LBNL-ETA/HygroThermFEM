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

    const std::vector<std::vector<double>> correctWaterContentSolution{{110.227073604, 110.227073604, 110.012394674, 110.012394674, 110.012394674, 110.012394674, 110.227073604, 110.227073604},
 {110.43212272, 110.43212272, 110.035369587, 110.035369587, 110.035369587, 110.035369587, 110.43212272, 110.43212272},
 {110.618344111, 110.61834411, 110.067366027, 110.067366027, 110.067366027, 110.067366027, 110.61834411, 110.61834411},
 {110.788480341, 110.788480341, 110.107046598, 110.107046598, 110.107046598, 110.107046598, 110.788480341, 110.788480341},
 {110.944881735, 110.944881735, 110.153265076, 110.153265076, 110.153265076, 110.153265076, 110.944881735, 110.944881735},
 {111.089560481, 111.089560481, 110.205040064, 110.205040064, 110.205040064, 110.205040064, 111.089560481, 111.089560481},
 {111.224237731, 111.224237731, 110.261532058, 110.261532058, 110.261532058, 110.261532058, 111.224237731, 111.224237731},
 {111.350396985, 111.350396986, 110.322017253, 110.322017253, 110.322017253, 110.322017253, 111.350396985, 111.350396985},
 {111.469277786, 111.469277786, 110.385891214, 110.385891214, 110.385891214, 110.385891214, 111.469277786, 111.469277786},
 {111.581953351, 111.581953351, 110.452630404, 110.452630404, 110.452630404, 110.452630404, 111.581953351, 111.581953351}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{1.99758110978e-12, 1.90217766996e-12, 1.82013380949e-12, 1.75034925809e-12, 1.69030928016e-12, 1.63859931893e-12, 1.59474782036e-12, 1.55702592475e-12, 1.52433344325e-12, 1.49651317909e-12};
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{78.8557392043, 78.8557392043, 65.4586114894, 65.4586114894, 65.4586114894, 65.4586114894, 78.8557392043, 78.8557392043},
 {85.8867106123, 85.8867106123, 66.6632540914, 66.6632540914, 66.6632540914, 66.6632540914, 85.8867106123, 85.8867106124},
 {89.2945311725, 89.2945311725, 68.3871255388, 68.3871255388, 68.3871255388, 68.3871255389, 89.2945311724, 89.2945311725},
 {90.9897409008, 90.9897409008, 70.32442862, 70.32442862, 70.32442862, 70.32442862, 90.9897409008, 90.9897409009},
 {91.9127735379, 91.9127735379, 72.2977789083, 72.2977789083, 72.2977789083, 72.2977789083, 91.9127735379, 91.9127735379},
 {92.4927618122, 92.4927618122, 74.2227026311, 74.2227026311, 74.2227026311, 74.2227026311, 92.4927618122, 92.4927618122},
 {92.9193056122, 92.9193056122, 76.0625905375, 76.0625905375, 76.0625905375, 76.0625905375, 92.9193056122, 92.9193056122},
 {93.2750103468, 93.2750103468, 77.8030041967, 77.8030041967, 77.8030041967, 77.8030041967, 93.2750103468, 93.2750103468},
 {93.5955111733, 93.5955111733, 79.4394763955, 79.4394763955, 79.4394763955, 79.4394763955, 93.5955111733, 93.5955111733},
 {93.8562745181, 93.8562745181, 81.0069281447, 81.0069281447, 81.0069281447, 81.0069281447, 93.8562745181, 93.8562745181}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{2.37656161023e-08, 1.35049391709e-08, 8.25582502982e-09, 5.73518992643e-09, 4.51137355062e-09, 3.87031662964e-09, 3.48616480615e-09, 3.21567449215e-09, 2.99683305365e-09, 2.79355842895e-09};
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{58.9855539536, 58.9855539536, 31.9496382205, 31.9496382205, 31.9496382205, 31.9496382205, 58.9855539536, 58.9855539536},
 {73.0272330997, 73.0272330997, 34.5256110684, 34.5256110684, 34.5256110684, 34.5256110684, 73.0272330997, 73.0272330997},
 {80.6354270535, 80.6354270535, 37.4330847932, 37.4330847932, 37.4330847932, 37.4330847932, 80.6354270535, 80.6354270535},
 {84.1318570155, 84.1318570155, 41.0704906234, 41.0704906234, 41.0704906234, 41.0704906234, 84.1318570155, 84.1318570155},
 {86.3661086423, 86.3661086423, 44.4916911717, 44.4916911717, 44.4916911717, 44.4916911717, 86.3661086423, 86.3661086423},
 {87.6515702862, 87.6515702862, 47.8997955163, 47.8997955163, 47.8997955163, 47.8997955163, 87.6515702862, 87.6515702862},
 {88.4645683517, 88.4645683517, 51.2358247163, 51.2358247163, 51.2358247163, 51.2358247163, 88.4645683517, 88.4645683517},
 {89.1002691705, 89.1002691705, 54.4194259443, 54.4194259443, 54.4194259443, 54.4194259443, 89.1002691705, 89.1002691705},
 {89.6647899601, 89.6647899601, 57.4245341687, 57.4245341687, 57.4245341687, 57.4245341687, 89.6647899601, 89.6647899601},
 {90.2038184514, 90.2038184513, 60.2379729593, 60.2379729593, 60.2379729593, 60.2379729593, 90.2038184514, 90.2038184513}};

    TestHelper::dumpGolden("correctWaterContentSolution", results.moisture.values);
    TestHelper::expectNear(correctWaterContentSolution, results.moisture.values, 1e-6);

    const std::vector<double> correctHumidityError{3.27345362488e-07, 6.69076106099e-08, 4.50374770073e-08, 4.35086495113e-08, 3.89524269714e-08, 1.48856543343e-08, 9.94226033534e-09, 9.22363196698e-09, 8.6312433323e-09, 8.08546802626e-09};
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

    const std::vector<std::vector<double>> correctWaterContentSolution{{355.000218234, 355.000218234, 355.000087304, 355.000087303, 355.000043657, 355.000043657, 355.000043657, 355.000043657, 355.000087304, 355.000087303, 355.000218234, 355.000218234},
 {355.000358372, 355.000358371, 355.000186319, 355.000186318, 355.000114648, 355.000114647, 355.000114648, 355.000114648, 355.000186318, 355.000186319, 355.000358372, 355.000358372},
 {355.000471418, 355.000471417, 355.000285169, 355.000285169, 355.000199347, 355.000199347, 355.000199347, 355.000199347, 355.000285169, 355.000285169, 355.000471418, 355.000471418},
 {355.000574309, 355.000574309, 355.000382761, 355.000382761, 355.00029038, 355.000290381, 355.000290381, 355.000290381, 355.000382761, 355.00038276, 355.00057431, 355.000574309},
 {355.00067318, 355.00067318, 355.00047954, 355.00047954, 355.000384235, 355.000384235, 355.000384235, 355.000384235, 355.00047954, 355.00047954, 355.00067318, 355.00067318},
 {355.000770403, 355.000770403, 355.000575908, 355.000575908, 355.000479323, 355.000479324, 355.000479323, 355.000479324, 355.000575908, 355.000575908, 355.000770404, 355.000770403},
 {355.000866936, 355.000866936, 355.000672085, 355.000672085, 355.000574946, 355.000574946, 355.000574946, 355.000574946, 355.000672085, 355.000672085, 355.000866937, 355.000866936},
 {355.000963177, 355.000963177, 355.000768176, 355.000768176, 355.000670799, 355.000670799, 355.000670799, 355.000670799, 355.000768176, 355.000768176, 355.000963178, 355.000963178},
 {355.001059293, 355.001059294, 355.00086423, 355.00086423, 355.00076675, 355.00076675, 355.000766751, 355.00076675, 355.00086423, 355.00086423, 355.001059294, 355.001059294},
 {355.001155356, 355.001155356, 355.000960267, 355.000960267, 355.000862744, 355.000862744, 355.000862744, 355.000862744, 355.000960267, 355.000960267, 355.001155357, 355.001155357},
 {355.001251396, 355.001251395, 355.001056296, 355.001056296, 355.000958755, 355.000958755, 355.000958755, 355.000958756, 355.001056296, 355.001056297, 355.001251396, 355.001251396},
 {355.001347425, 355.001347425, 355.001152321, 355.001152322, 355.001054773, 355.001054774, 355.001054774, 355.001054774, 355.001152322, 355.001152322, 355.001347425, 355.001347425},
 {355.001443449, 355.001443449, 355.001248345, 355.001248345, 355.001150794, 355.001150795, 355.001150795, 355.001150795, 355.001248346, 355.001248346, 355.001443449, 355.001443449},
 {355.00153947, 355.001539471, 355.001344368, 355.001344368, 355.001246817, 355.001246817, 355.001246817, 355.001246817, 355.001344368, 355.001344368, 355.001539471, 355.001539471},
 {355.001635491, 355.001635491, 355.00144039, 355.00144039, 355.001342839, 355.001342839, 355.001342839, 355.001342839, 355.00144039, 355.00144039, 355.001635491, 355.001635491},
 {355.00173151, 355.00173151, 355.001536411, 355.001536411, 355.001438861, 355.001438861, 355.001438861, 355.001438861, 355.001536411, 355.001536411, 355.00173151, 355.00173151},
 {355.001827529, 355.001827529, 355.001632432, 355.001632431, 355.001534882, 355.001534882, 355.001534882, 355.001534882, 355.001632431, 355.001632431, 355.001827529, 355.001827528},
 {355.001923547, 355.001923547, 355.001728451, 355.001728451, 355.001630902, 355.001630902, 355.001630903, 355.001630902, 355.001728451, 355.001728451, 355.001923547, 355.001923546},
 {355.002019563, 355.002019564, 355.00182447, 355.00182447, 355.001726922, 355.001726922, 355.001726922, 355.001726922, 355.00182447, 355.00182447, 355.002019564, 355.002019564},
 {355.00211558, 355.00211558, 355.001920488, 355.001920488, 355.001822942, 355.001822941, 355.001822942, 355.001822941, 355.001920488, 355.001920489, 355.00211558, 355.002115581}};

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
