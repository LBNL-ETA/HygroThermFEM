#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

TEST(TestSaturationFunctions, BoundarySaturation)
{
    constexpr auto temperature = 20.0;

    const auto result = HygroThermFEM::vaporPressureAtTemperature(temperature);

    constexpr auto correctResult = 2331.215017;

    EXPECT_NEAR(correctResult, result, 1e-6);
}

TEST(TestSaturationFunctions, Saturation)
{
    constexpr auto temperature = 20.0;

	const auto result = HygroThermFEM::saturationConcentrationAtTemperature(temperature);

	constexpr auto correctResult = 0.017235141;

	EXPECT_NEAR(correctResult, result, 1e-8);
}