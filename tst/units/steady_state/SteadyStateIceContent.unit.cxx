#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

#include "BeamBuilder.hxx"
#include "TestMaterials.hxx"

// The steady solve keeps the latent-heat capacity out of the equations, but the
// reported liquid/ice split must still follow the final temperature: condensed water
// below the freezing ramp is ice, above it liquid. Without the split the steady ice
// field would be identically zero however cold the model.
TEST(SteadyStateIceContent, ColdFaceReportsIce)
{
    HygroThermFEM::MultiDomain multiDomain({.performThermal = true, .performMoisture = false});

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

    // Steady faces at about -13 C and +11 C: the split must flip between the two ends.
    const HygroThermFEM::FixedBCHCCoefficients coldSurface{-18.0, 20.0, 0.0};
    const HygroThermFEM::FixedBCHCCoefficients warmSurface{21.0, 10.0, 0.0};
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Left, coldSurface);
    builder.applyBC_FixedHc(TestHelper::BeamBuilder::Edge::Right, warmSurface);

    const auto solution = multiDomain.steadyState();
    ASSERT_TRUE(solution.converged);
    ASSERT_EQ(solution.iceContent.size(), solution.temperature.size());

    // Column-major nodes, one element row: node 0 sits on the cold face, the last node
    // on the warm face.
    const std::size_t coldFace{0u};
    const std::size_t warmFace{solution.temperature.size() - 1u};
    ASSERT_LT(solution.temperature[coldFace], -0.1);
    ASSERT_GT(solution.temperature[warmFace], 0.0);

    const double condensedAtColdFace =
      solution.waterContent[coldFace] - solution.vaporContent[coldFace];
    EXPECT_GT(condensedAtColdFace, 50.0);
    EXPECT_NEAR(solution.iceContent[coldFace], condensedAtColdFace, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[coldFace], 0.0, 1e-6);

    const double condensedAtWarmFace =
      solution.waterContent[warmFace] - solution.vaporContent[warmFace];
    EXPECT_NEAR(solution.iceContent[warmFace], 0.0, 1e-6);
    EXPECT_NEAR(solution.liquidWaterContent[warmFace], condensedAtWarmFace, 1e-6);
}
