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
	const FenestrationCommon::TabularFunction curve{ { 1, 10 },
																				 { 2, 20 },
																				 { 3, 30 } };

	double interpolationPoint = 2.5;

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 25, result, 1e-6 );

}

TEST_F( CurveTest, TestFirstDerivative ) {
	SCOPED_TRACE( "Begin Test: Test first derivative of tabular linear curve." );
	const FenestrationCommon::FirstDerivativeCurve curve{ { 1, 10 },
																												{ 2, 20 },
																												{ 3, 30 } };

	double interpolationPoint = 2.5;

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 10, result, 1e-6 );
}

TEST_F( CurveTest, TestTabularLogarithmic ) {
	SCOPED_TRACE( "Begin Test: Test tabular logarithmic curve." );
	const FenestrationCommon::TabularFunction curve( { { 1, 10 },
																					 { 2, 20 },
																					 { 3, 30 } },
																				 FenestrationCommon::Interpolation::Logarithmic );

	double interpolationPoint = 2.5;

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 24.4948974, result, 1e-6 );

}

TEST_F( CurveTest, TestSuctionCurve ) {
	SCOPED_TRACE( "Begin Test: Test suction curve." );
	const FenestrationCommon::SuctionCurve curve( { { 1, 10 },
																									{ 2, 20 },
																									{ 3, 30 } } );

	/// First segment should have constant values
	double interpolationPoint = 1.5;
	auto result = curve.value( interpolationPoint );
	EXPECT_NEAR( 10, result, 1e-6 );

	/// Test outside of curve
	interpolationPoint = 0.5;
	result = curve.value( interpolationPoint );
	EXPECT_NEAR( 10, result, 1e-6 );

	/// Other segments should have logarithmic interpolation
	interpolationPoint = 2.5;
	result = curve.value( interpolationPoint );
	EXPECT_NEAR( 24.4948974, result, 1e-6 );

}

TEST_F( CurveTest, TestConstantCurve ) {
	SCOPED_TRACE( "Begin Test: Test tabular logarithmic curve." );
	const FenestrationCommon::Constant cons( 5 );

	double interpolationPoint = 2.5;
	auto result = cons.value( interpolationPoint );
	EXPECT_NEAR( 5, result, 1e-6 );

}

TEST_F( CurveTest, TestTabularOutOfRangeBack ) {
	SCOPED_TRACE( "Begin Test: Test tabular out of range." );
	const FenestrationCommon::TabularFunction curve{ { 1, 10 },
																				 { 2, 20 },
																				 { 3, 30 } };

	double interpolationPoint = 3.5;
	auto result = curve.value( interpolationPoint );
	EXPECT_NEAR( 30, result, 1e-6 );

}

TEST_F( CurveTest, TestTabularOutOfRangeFront ) {
	SCOPED_TRACE( "Begin Test: Test tabular out of range." );
	const FenestrationCommon::TabularFunction curve{ { 1, 10 },
																				 { 2, 20 },
																				 { 3, 30 } };

	double interpolationPoint = 0.5;
	auto result = curve.value( interpolationPoint );
	EXPECT_NEAR( 10, result, 1e-6 );

}
