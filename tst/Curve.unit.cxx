#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

class CurveTest : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {

	}

};

TEST_F( CurveTest, TestTabularLinear ) {
	SCOPED_TRACE( "Begin Test: Test tabular linear curve." );
	const FenestrationCommon::Curve curve{ { 1, 10 },
																				 { 2, 20 },
																				 { 3, 30 } };

	double interpolationPoint = 2.5;

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 25, result, 1e-6 );

	result = curve.firstDerivative( interpolationPoint );

	EXPECT_NEAR( 10, result, 1e-6 );

}

TEST_F( CurveTest, TestTabularLogarithmic ) {
	SCOPED_TRACE( "Begin Test: Test tabular logarithmic curve." );
	const FenestrationCommon::Curve curve( { { 1, 10 },
																					 { 2, 20 },
																					 { 3, 30 } },
																				 FenestrationCommon::Interpolator::Logarithmic );

	double interpolationPoint = 2.5;

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 24.4948974, result, 1e-6 );

	result = curve.firstDerivative( interpolationPoint );

	EXPECT_NEAR( 9.931826, result, 1e-6 );

}

TEST_F( CurveTest, TestConstantCurve ) {
	SCOPED_TRACE( "Begin Test: Test tabular logarithmic curve." );
	const FenestrationCommon::Constant cons( 5 );

	double interpolationPoint = 2.5;

	auto result = cons.value( interpolationPoint );

	EXPECT_NEAR( 5, result, 1e-6 );

	result = cons.firstDerivative( interpolationPoint );

	EXPECT_NEAR( 0, result, 1e-6 );

}
