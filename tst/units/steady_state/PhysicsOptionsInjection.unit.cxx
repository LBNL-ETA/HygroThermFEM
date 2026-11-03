#include <cmath>
#include <optional>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::PhysicsOptions;
using HygroThermFEM::SimulationProperties;
using HygroThermFEM::State;

namespace
{
    void createNodes(HygroThermFEM::MultiDomain & multiDomain)
    {
        constexpr double initialTemperature = 21.0;
        constexpr double initialPressure = 101325.0;
        constexpr auto liquidPercent = 1.0;

        const auto makeState = [&](const double humidity) {
            return State{.temperature = initialTemperature,
                         .humidity = humidity,
                         .pressure = initialPressure,
                         .liquidPercent = liquidPercent};
        };

        multiDomain.nodes().createNode({.index = 1, .x = 1, .y = 5, .state = makeState(0.0)});
        multiDomain.nodes().createNode({.index = 2, .x = 1, .y = 0, .state = makeState(0.0)});
        multiDomain.nodes().createNode({.index = 3, .x = 0.5, .y = 5, .state = makeState(0.5)});
        multiDomain.nodes().createNode({.index = 4, .x = 0.5, .y = 0, .state = makeState(0.5)});
        multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5, .state = makeState(1.0)});
        multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = makeState(1.0)});
    }

    //! Solves the SteadyState_2D_1-style two-element model, optionally with injected physics.
    //! Without injection the domain falls back to the global SimulationProperties singleton.
    HygroThermFEM::Solution solveModel(const std::optional<PhysicsOptions> & physics)
    {
        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});
        if(physics.has_value())
        {
            multiDomain.setPhysicsOptions(physics.value());
        }

        createNodes(multiDomain);

        auto params = TestHelper::TestMaterial();
        params.porosity = 0.18;
        const auto & material = multiDomain.materials().createSolidMaterial(params);

        multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
        multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

        constexpr auto hcPin = 1e20;
        const HygroThermFEM::FixedBCHCCoefficients bcCoeff1{0.0, hcPin, 0.8};
        const HygroThermFEM::FixedBCHCCoefficients bcCoeff2{20.0, hcPin, 0.0};
        multiDomain.createBC_FixedHc(1, 2, bcCoeff1);
        multiDomain.createBC_FixedHc(6, 5, bcCoeff2);

        return multiDomain.steadyState();
    }

    void expectSameSolution(const HygroThermFEM::Solution & lhs,
                            const HygroThermFEM::Solution & rhs)
    {
        ASSERT_EQ(lhs.temperature.size(), rhs.temperature.size());
        for(auto idx = 0u; idx < lhs.temperature.size(); ++idx)
        {
            EXPECT_NEAR(lhs.temperature[idx], rhs.temperature[idx], 1e-9);
        }
        ASSERT_EQ(lhs.humidity.size(), rhs.humidity.size());
        for(auto idx = 0u; idx < lhs.humidity.size(); ++idx)
        {
            EXPECT_NEAR(lhs.humidity[idx], rhs.humidity[idx], 1e-9);
        }
    }
}   // namespace

class PhysicsOptionsInjection : public testing::Test
{
protected:
    void TearDown() override
    {
        SimulationProperties::Instance().reset();
    }
};

// The injected options must win over the global singleton: with the globals dirtied to exclude
// heat of evaporation, a domain injected with the engine defaults must reproduce the
// clean-globals solution, not the dirty-globals one.
TEST_F(PhysicsOptionsInjection, InjectedDefaultsIgnoreDirtyGlobals)
{
    SimulationProperties::Instance().setCalculationParameters(
      /*excludeWaterLiquidTransportation=*/false,
      /*excludeHeatOfEvaporation=*/true,
      /*excludeCapillaryConduction=*/false,
      /*excludeVaporDiffusionConduction=*/false,
      /*thermalConductivityMoistureAndTemperatureDependent=*/false);

    const auto injected = solveModel(PhysicsOptions{});
    const auto dirtyGlobals = solveModel(std::nullopt);

    SimulationProperties::Instance().reset();
    const auto cleanGlobals = solveModel(std::nullopt);

    expectSameSolution(injected, cleanGlobals);

    // The dirty globals would have changed the interior temperature, so an injected domain
    // that silently read the singleton would fail the comparison above. The solver is
    // deterministic, so any threshold above its 1e-5 iteration tolerance is a physics
    // difference, not numerical noise.
    EXPECT_GT(std::abs(injected.temperature[2] - dirtyGlobals.temperature[2]), 1e-4);
}

// Two domains in the same process solve with DIFFERENT physics: the per-instance options are
// what makes concurrent solves with different physics safe.
TEST_F(PhysicsOptionsInjection, TwoDomainsWithDifferentInjectedPhysics)
{
    const auto evaporationExcluded = solveModel(PhysicsOptions{.excludeHeatOfEvaporation = true});
    const auto engineDefaults = solveModel(PhysicsOptions{});

    ASSERT_EQ(evaporationExcluded.temperature.size(), engineDefaults.temperature.size());
    EXPECT_GT(std::abs(evaporationExcluded.temperature[2] - engineDefaults.temperature[2]), 1e-4);
}
