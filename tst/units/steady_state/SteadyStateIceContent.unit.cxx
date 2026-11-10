#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

namespace
{
    //! A wet stucco beam whose steady faces settle at about -13 C and +11 C, solved with
    //! the given physics. Column-major nodes, one element row: node 0 sits on the cold
    //! face, the last node on the warm face.
    HygroThermFEM::Solution solveColdToWarmBeam(const HygroThermFEM::PhysicsOptions & physics)
    {
        HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});
        multiDomain.setPhysicsOptions(physics);

        const auto & stucco = multiDomain.materials().createSolidMaterial(TestHelper::Stucco());

        const HygroThermFEM::State initialState({.temperature = 10.0,
                                                 .humidity = 0.99,
                                                 .pressure = 101325.0,
                                                 .liquidPercent = 1.0});

        TestHelper::BeamBuilder builder(multiDomain);
        builder.xStart(0.0)
          .height(0.05)
          .numElementsY(1)
          .state(initialState)
          .addSegment({.material = stucco.name(), .numElementsX = 20, .width = 0.2})
          .build();

        const HygroThermFEM::FixedBCHCCoefficients coldSurface{-18.0, 20.0, 0.0};
        const HygroThermFEM::FixedBCHCCoefficients warmSurface{21.0, 10.0, 0.0};
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, coldSurface);
        builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, warmSurface);

        return multiDomain.steadyState();
    }

    double condensedWater(const HygroThermFEM::Solution & solution, const std::size_t node)
    {
        return solution.waterContent[node] - solution.vaporContent[node];
    }
}   // namespace

// The steady solve keeps the latent-heat capacity out of the equations, but with freezing
// modelled the reported liquid/ice split must still follow the final temperature:
// condensed water below the freezing ramp is ice, above it liquid.
TEST(SteadyStateIceContent, ColdFaceReportsIceWhenFreezingIsModelled)
{
    const auto solution = solveColdToWarmBeam({.excludeLatentHeatOfFusion = false});
    ASSERT_TRUE(solution.converged);
    ASSERT_EQ(solution.iceContent.size(), solution.temperature.size());

    const std::size_t coldFace{0u};
    const std::size_t warmFace{solution.temperature.size() - 1u};
    ASSERT_LT(solution.temperature[coldFace], -0.1);
    ASSERT_GT(solution.temperature[warmFace], 0.0);

    const double condensedAtColdFace{condensedWater(solution, coldFace)};
    EXPECT_GT(condensedAtColdFace, 50.0);
    EXPECT_NEAR(solution.iceContent[coldFace], condensedAtColdFace, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[coldFace], 0.0, 1e-6);

    EXPECT_NEAR(solution.iceContent[warmFace], 0.0, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[warmFace], condensedWater(solution, warmFace), 1e-6);
}

// With the latent heat of fusion excluded freezing is not modelled, in a steady solve as in
// a transient one: all condensed water is reported as liquid however cold the model.
TEST(SteadyStateIceContent, NoIceWhenFreezingIsExcluded)
{
    const auto solution = solveColdToWarmBeam({.excludeLatentHeatOfFusion = true});
    ASSERT_TRUE(solution.converged);

    const std::size_t coldFace{0u};
    ASSERT_LT(solution.temperature[coldFace], -0.1);

    const double condensedAtColdFace{condensedWater(solution, coldFace)};
    EXPECT_GT(condensedAtColdFace, 50.0);
    EXPECT_NEAR(solution.iceContent[coldFace], 0.0, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[coldFace], condensedAtColdFace, 1e-6);
}

// The option changes only what is reported: the steady temperatures are the same.
TEST(SteadyStateIceContent, TemperaturesDoNotDependOnTheOption)
{
    const auto modelled = solveColdToWarmBeam({.excludeLatentHeatOfFusion = false});
    const auto excluded = solveColdToWarmBeam({.excludeLatentHeatOfFusion = true});

    ASSERT_EQ(modelled.temperature.size(), excluded.temperature.size());
    for(std::size_t node = 0u; node < modelled.temperature.size(); ++node)
    {
        EXPECT_NEAR(modelled.temperature[node], excluded.temperature[node], 1e-6);
    }
}
