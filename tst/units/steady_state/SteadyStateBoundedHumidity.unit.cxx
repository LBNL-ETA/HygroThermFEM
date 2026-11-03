#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"
#include "TestMaterials.hxx"
#include "TestHelpers.hxx"

using HygroThermFEM::Nodes;
using HygroThermFEM::State;

// A steady vapor-diffusion solve across a temperature gradient supersaturates where the local
// saturation concentration drops faster than the vapor concentration (Glaser): with the warm
// face at 20 C / RH 1.0 and the cold face at 0 C / RH 0.99, the unconstrained mid-plane
// humidity lands around 1.09. The bounded solve must pin such nodes at saturation (RH 1.0)
// instead -- the same ceiling the transient solver enforces through its clamped Newton
// corrections -- so no node may ever report humidity above 1.
TEST(SteadyStateBoundedHumidity, CondensationZonePinsAtSaturation)
{
    constexpr double initialTemperature = 10.0;
    constexpr double initialHumidity = 0.5;
    constexpr double initialPressure = 101325.0;
    constexpr auto liquidPercent = 1.0;

    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = true});

    const State initialState{.temperature = initialTemperature,
                             .humidity = initialHumidity,
                             .pressure = initialPressure,
                             .liquidPercent = liquidPercent};

    multiDomain.nodes().createNode({.index = 1, .x = 1, .y = 5, .state = initialState});
    multiDomain.nodes().createNode({.index = 2, .x = 1, .y = 0, .state = initialState});
    multiDomain.nodes().createNode({.index = 3, .x = 0.5, .y = 5, .state = initialState});
    multiDomain.nodes().createNode({.index = 4, .x = 0.5, .y = 0, .state = initialState});
    multiDomain.nodes().createNode({.index = 5, .x = 0, .y = 5, .state = initialState});
    multiDomain.nodes().createNode({.index = 6, .x = 0, .y = 0, .state = initialState});

    auto params = TestHelper::TestMaterial();
    params.porosity = 0.18;
    // Zero-valued liquid transport keeps the moisture equation vapor-only, so the solve is a
    // single deterministic linear system and the supersaturation cannot be redistributed away --
    // the projection itself is what gets exercised.
    params.liquidTransportCurve = {{0, 0}, {180, 0}};
    params.sorptionCurve = {{0, 0}, {1, 180}};
    const auto & material = multiDomain.materials().createSolidMaterial(params);

    multiDomain.createElement({.node1 = 3, .node2 = 4, .node3 = 2, .node4 = 1, .material = material.name()});
    multiDomain.createElement({.node1 = 6, .node2 = 4, .node3 = 3, .node4 = 5, .material = material.name()});

    // Dirichlet-like faces through huge film coefficients: cold face (nodes 1-2) at 0 C and
    // RH 0.99; warm face (nodes 5-6) at 20 C and RH 1.0.
    constexpr auto boundaryFilmCoefficient = 1e20;
    const HygroThermFEM::FixedBCHCCoefficients coldFace{0.0, boundaryFilmCoefficient, 0.99};
    const HygroThermFEM::FixedBCHCCoefficients warmFace{20.0, boundaryFilmCoefficient, 1.0};

    multiDomain.createBC_FixedHc(1, 2, coldFace);
    multiDomain.createBC_FixedHc(6, 5, warmFace);

    const auto solution = multiDomain.steadyState();
    const auto & humidity = solution.humidity;

    TestHelper::dumpGolden("boundedHumidity", humidity);
    ASSERT_EQ(humidity.size(), 6u);

    // No node may exceed saturation.
    constexpr double boundTolerance = 1e-9;
    for(const auto value : humidity)
    {
        EXPECT_LE(value, 1.0 + boundTolerance);
        EXPECT_GE(value, 0.0 - boundTolerance);
    }

    // The faces hold their boundary humidities...
    EXPECT_NEAR(humidity[0], 0.99, 1e-6);
    EXPECT_NEAR(humidity[1], 0.99, 1e-6);
    EXPECT_NEAR(humidity[4], 1.0, 1e-6);
    EXPECT_NEAR(humidity[5], 1.0, 1e-6);

    // ...and the mid-plane, which the vapor-only unconstrained solve drives above saturation,
    // sits pinned at exactly 1 -- without pinned interior nodes no supersaturation occurred and
    // the test would no longer exercise the projection.
    EXPECT_NEAR(humidity[2], 1.0, 1e-6);
    EXPECT_NEAR(humidity[3], 1.0, 1e-6);
}
