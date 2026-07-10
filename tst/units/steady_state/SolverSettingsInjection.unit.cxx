#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::State;

// Exercises the injected SolverSettings path (MultiDomain::setSolverSettings) end-to-end: the same
// model as SteadyState_2D_1, but the solver configuration is supplied explicitly instead of being
// read from the global SimulationProperties / Timesteps singletons. Injecting the engine-default
// values must reproduce the reference solution, proving the injected configuration is threaded
// through the solve rather than ignored.
TEST(SolverSettingsInjection, SteadyStateWithInjectedDefaults)
{
    constexpr double initialTemperature = 21.0;
    constexpr double initialPressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

    // Engine-default solver configuration, supplied explicitly rather than via the singletons.
    multiDomain.setSolverSettings({.relaxationParameter = 1.0,
                                   .errorTolerance = 1e-5,
                                   .maxNumberOfIterations = 50,
                                   .maxDivisions = 3,
                                   .numberOfSubtimesteps = 10});

    multiDomain.nodes().createNode({.index = 1, .x = 1, .y = 5, .state = State{
        .temperature = initialTemperature, .humidity = 0, .pressure = initialPressure, .liquidPercent = liquidPercent}});
    multiDomain.nodes().createNode({.index = 2, .x = 1, .y = 0, .state = State{
        .temperature = initialTemperature, .humidity = 0, .pressure = initialPressure, .liquidPercent = liquidPercent}});
    multiDomain.nodes().createNode({.index = 3, .x = 0.5, .y = 5, .state = State{
        .temperature = initialTemperature, .humidity = 0.5, .pressure = initialPressure, .liquidPercent = liquidPercent}});
    multiDomain.nodes().createNode({.index = 4, .x = 0.5, .y = 0, .state = State{
        .temperature = initialTemperature, .humidity = 0.5, .pressure = initialPressure, .liquidPercent = liquidPercent}});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5, .state = State{
        .temperature = initialTemperature, .humidity = 1, .pressure = initialPressure, .liquidPercent = liquidPercent}});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = State{
        .temperature = initialTemperature, .humidity = 1, .pressure = initialPressure, .liquidPercent = liquidPercent}});

    auto params = TestHelper::TestMaterial();
    params.porosity = 0.18;
    params.liquidTransportCurve = {{0, 0}, {180, 7E-7}};
    params.sorptionCurve = {{0, 0}, {1, 180}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    constexpr auto hc = 1e20;
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{0.0, hc, 0.0};
    const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{20.0, hc, 1.0};
    multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
    multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

    const auto solution = multiDomain.steadyState();

    const std::vector<double> correctTemperature{2.13180767e-19, 2.13180767e-19, 10.6590052, 10.6590052, 20, 20};
    TestHelper::dumpGolden("correctTemperature", solution.temperature);
    ASSERT_EQ(solution.temperature.size(), correctTemperature.size());
    for(auto i = 0u; i < correctTemperature.size(); ++i)
    {
        EXPECT_NEAR(solution.temperature[i], correctTemperature[i], 1e-6);
    }
}
