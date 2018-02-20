#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using MoisThermFEM::Property;
using MoisThermFEM::State;

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
	SCOPED_TRACE( "Begin Test: Test tabular linear." );
	const MoisThermFEM::TabularFunction curve( { { 1, 10 },
																							 { 2, 20 },
																							 { 3, 30 } }, Property::temperature );

	State interpolationPoint( 2.5, 0, 101325 );

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 25, result, 1e-6 );

	auto max = curve.max();
	EXPECT_NEAR( 30, max, 1e-6 );

	auto min = curve.min();
	EXPECT_NEAR( 10, min, 1e-6 );

}

TEST_F( CurveTest, TestFirstDerivative ) {
	SCOPED_TRACE( "Begin Test: Test first derivative of tabular linear curve." );
	const MoisThermFEM::FirstDerivativeFunction curve( { { 1, 10 },
																											 { 2, 20 },
																											 { 3, 30 } }, Property::temperature );

	State interpolationPoint( 2.5, 0, 101325 );

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 10, result, 1e-6 );
}

TEST_F( CurveTest, TestTabularLogarithmic ) {
	SCOPED_TRACE( "Begin Test: Test tabular logarithmic curve." );
	const MoisThermFEM::TabularFunction curve( { { 1, 10 },
																							 { 2, 20 },
																							 { 3, 30 } },
																						 Property::temperature,
																						 FenestrationCommon::Interpolation::Logarithmic );

	State interpolationPoint( 2.5, 0, 101325 );

	auto result = curve.value( interpolationPoint );

	EXPECT_NEAR( 24.4948974, result, 1e-6 );

}

TEST_F( CurveTest, TestSuctionCurve ) {
	SCOPED_TRACE( "Begin Test: Test suction function." );
	const MoisThermFEM::SuctionFunction curve( { { 1, 10 },
																							 { 2, 20 },
																							 { 3, 30 } }, Property::temperature );

	/// First segment should have constant values
	State interpolationPoint( 1.5, 0, 101325 );
	auto result = curve.value( interpolationPoint );
	EXPECT_NEAR( 10, result, 1e-6 );

	/// Test outside of curve
	State interpolationPoint1( 0.5, 0, 101325 );
	result = curve.value( interpolationPoint1 );
	EXPECT_NEAR( 10, result, 1e-6 );

	/// Other segments should have logarithmic interpolation
	State interpolationPoint2( 2.5, 0, 101325 );
	result = curve.value( interpolationPoint2 );
	EXPECT_NEAR( 24.4948974, result, 1e-6 );

}

TEST_F( CurveTest, TestConstantCurve ) {
	SCOPED_TRACE( "Begin Test: Test tabular logarithmic." );
	const MoisThermFEM::Constant cons( 5 );

	State interpolationPoint( 2.5, 0, 101325 );
	auto result = cons.value( interpolationPoint );
	EXPECT_NEAR( 5, result, 1e-6 );

}

TEST_F( CurveTest, TestTabularOutOfRangeBack ) {
	SCOPED_TRACE( "Begin Test: Test tabular out of range." );
	const MoisThermFEM::TabularFunction curve( { { 1, 10 },
																							 { 2, 20 },
																							 { 3, 30 } }, Property::temperature );

	State interpolationPoint( 3.5, 0, 101325 );
	auto result = curve.value( interpolationPoint );
	EXPECT_NEAR( 30, result, 1e-6 );

}

TEST_F( CurveTest, TestTabularOutOfRangeFront ) {
	SCOPED_TRACE( "Begin Test: Test tabular out of range." );
	const MoisThermFEM::TabularFunction curve( { { 1, 10 },
																							 { 2, 20 },
																							 { 3, 30 } }, Property::temperature );

	State interpolationPoint( 0.5, 0, 101325 );
	auto result = curve.value( interpolationPoint );
	EXPECT_NEAR( 10, result, 1e-6 );

}

TEST_F( CurveTest, TestComposition1 ) {
	SCOPED_TRACE( "Begin Test: Composition (multiplication) of two functions." );
	std::unique_ptr< MoisThermFEM::IFunction > cons = fem::make_unique< MoisThermFEM::Constant >( 5 );
	const MoisThermFEM::TabularFunction tabular( { { 1, 10 },
																								 { 2, 20 },
																								 { 3, 30 } }, Property::temperature, cons );

	State interpolationPoint( 2.5, 0, 101325 );

	auto result = tabular.value( interpolationPoint );

	EXPECT_NEAR( 125, result, 1e-6 );

}

TEST_F( CurveTest, TestPorosityCalculation ) {
	SCOPED_TRACE( "Begin Test: Calculate liquid and air porosities." );

	std::unique_ptr< MoisThermFEM::IFunction > waterContent =
			std::unique_ptr< MoisThermFEM::IFunction >(
					new MoisThermFEM::TabularFunction( { { 0.000, 0.0 },
																							 { 0.500, 0.5 },
																							 { 0.800, 1.4 },
																							 { 0.900, 2.6 },
																							 { 0.930, 3.6 },
																							 { 0.950, 4.7 },
																							 { 0.970, 7.1 },
																							 { 0.990, 14.8 },
																							 { 0.995, 20.9 },
																							 { 0.999, 33.0 },
																							 { 1.000, 40.0 } }, Property::humidity ) );


	auto maxWaterContent = waterContent->value( State( 0, 1, 0 ) );
	const auto materialPorosity = 0.05;

	std::unique_ptr< MoisThermFEM::IFunction > waterFill = fem::make_unique< MoisThermFEM::Constant >(
			materialPorosity / maxWaterContent, waterContent );

	State outdoor( 10, 0.98, 101325 );

	auto result = waterFill->value( outdoor );
	EXPECT_NEAR( 0.0136875, result, 1e-6 );

	const MoisThermFEM::Constant airFill( materialPorosity, waterFill, MoisThermFEM::Operation::SUB );

	result = airFill.value( outdoor );
	EXPECT_NEAR( 0.0363125, result, 1e-6 );

}

TEST_F( CurveTest, TestSaturationFunction ) {
	SCOPED_TRACE( "Begin Test: Test saturation function." );
	const MoisThermFEM::SaturationFunction sat( Property::temperature );

	State interpolationPoint( 293.15, 0, 101325 );
	auto result = sat.value( interpolationPoint );
	EXPECT_NEAR( 0.017235141, result, 1e-6 );

}
