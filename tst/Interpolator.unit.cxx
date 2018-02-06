#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using FenestrationCommon::InterpolatorFactory;
using FenestrationCommon::Interpolator;

class InterpolatorTest : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {

	}

};

TEST_F( InterpolatorTest, TestLinear ) {
	SCOPED_TRACE( "Begin Test: Linear integrator." );
	const auto interpolator = InterpolatorFactory::getInterpolator( Interpolator::Linear );

	std::pair< double, double > pt1{ 31, 4 };
	std::pair< double, double > pt2{ 52, 8 };
	double interpolationPoint = 40;


	/// Now test commutation
	auto result = interpolator->interpolate( pt1, pt2, interpolationPoint );

	EXPECT_NEAR( 5.714285714, result, 1e-6 );

}

TEST_F( InterpolatorTest, TestLogarithmic ) {
	SCOPED_TRACE( "Begin Test: Logarithmic integrator." );

	const auto interpolator = InterpolatorFactory::getInterpolator( Interpolator::Logarithmic );

	std::pair< double, double > pt1{ 31, 4 };
	std::pair< double, double > pt2{ 52, 8 };
	double interpolationPoint = 40;


	/// Now test commutation
	auto result = interpolator->interpolate( pt1, pt2, interpolationPoint );

	EXPECT_NEAR( 5.383600771, result, 1e-6 );

}
