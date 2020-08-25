#include <memory>
#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

class TestSaturationFunctions : public testing::Test
{
protected:
    void SetUp() override
    {}
};

TEST_F(TestSaturationFunctions, BoundarySaturation)
{
    const auto temperature = 20.0;

    const auto result = HygroThermFEM::vaporPressureAtTemperature(temperature);

    const auto correctResult = 2331.215017;

    EXPECT_NEAR(correctResult, result, 1e-6);
}

TEST_F(TestSaturationFunctions, Saturation)
{
    const auto temperature = 20.0;

	const auto result = HygroThermFEM::saturationConcentrationAtTemperature(temperature);

	const auto correctResult = 0.017235141;

	EXPECT_NEAR(correctResult, result, 1e-8);
}