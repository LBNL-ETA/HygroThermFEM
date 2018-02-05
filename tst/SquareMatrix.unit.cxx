#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using FenestrationCommon::SquareMatrix;

/////////////////////////////////////////////////////////////////////////////////////
/// Transient heat transfer example on Sandstone specimen using data from database
///   Lumped mass matrix
///   Time-step 1 hour
///   Six nodes block at initial temperatures in nodes of 100 degrees
///   Initial temperature boundary conditions at nodes 5 and 6 are 12 degrees
///   Solution achieved with linear solver (no iterations required in this case
/////////////////////////////////////////////////////////////////////////////////////

class SquareMatrixTest : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {

	}

};

TEST_F( SquareMatrixTest, TestAddition ) {
	SCOPED_TRACE( "Begin Test: Matrix addition." );
	SquareMatrix< double > a{ { 1, 2 },
														{ 3, 4 } };

	SquareMatrix< double > b{ { 2, 8 },
														{ 3, 5 } };

	auto result = a + b;

	SquareMatrix< double > correctResult{ { 3, 10 },
																				{ 6, 9 } };

	EXPECT_EQ( correctResult.size(), result.size() );

	for ( auto i = 0u; i < correctResult.size(); ++i ) {
		for ( auto j = 0u; j < correctResult.size(); ++j ) {
			EXPECT_NEAR( correctResult[ i ][ j ], result[ i ][ j ], 1e-6 );
		}
	}

}

TEST_F( SquareMatrixTest, TestMultiplication ) {
	SCOPED_TRACE( "Begin Test: Matrix multiplication." );
	SquareMatrix< double > a{ { 1, 2 },
														{ 3, 4 } };

	SquareMatrix< double > b{ { 2, 8 },
														{ 3, 5 } };

	auto result = a * b;

	SquareMatrix< double > correctResult{ { 8, 18 },
																				{ 18, 44 } };

	EXPECT_EQ( correctResult.size(), result.size() );

	for ( auto i = 0u; i < correctResult.size(); ++i ) {
		for ( auto j = 0u; j < correctResult.size(); ++j ) {
			EXPECT_NEAR( correctResult[ i ][ j ], result[ i ][ j ], 1e-6 );
		}
	}

	/// Now test inverse multiplication
	result = b * a;

	SquareMatrix< double > correctInverse{ { 26, 36 },
																				 { 18, 26 } };

	EXPECT_EQ( correctInverse.size(), result.size() );

	for ( auto i = 0u; i < correctInverse.size(); ++i ) {
		for ( auto j = 0u; j < correctInverse.size(); ++j ) {
			EXPECT_NEAR( correctInverse[ i ][ j ], result[ i ][ j ], 1e-6 );
		}
	}
}

TEST_F( SquareMatrixTest, TestMultiplicationBySingleValue ) {
	SCOPED_TRACE( "Begin Test: Matrix multiplication by single value." );
	SquareMatrix< double > a{ { 1, 2 },
														{ 3, 4 } };

	double b = 2;

	auto result = a * b;

	SquareMatrix< double > correctResult{ { 2, 4 },
																				{ 6, 8 } };

	EXPECT_EQ( correctResult.size(), result.size() );

	for ( auto i = 0u; i < correctResult.size(); ++i ) {
		for ( auto j = 0u; j < correctResult.size(); ++j ) {
			EXPECT_NEAR( correctResult[ i ][ j ], result[ i ][ j ], 1e-6 );
		}
	}

	/// Now test inverse multiplication
	result = b * a;

	EXPECT_EQ( correctResult.size(), result.size() );

	for ( auto i = 0u; i < correctResult.size(); ++i ) {
		for ( auto j = 0u; j < correctResult.size(); ++j ) {
			EXPECT_NEAR( correctResult[ i ][ j ], result[ i ][ j ], 1e-6 );
		}
	}
}

TEST_F( SquareMatrixTest, TestMultiplicationWithVectors ) {
	SCOPED_TRACE( "Begin Test: Matrix multiplication with vector." );
	SquareMatrix< double > a{ { 1, 2 },
														{ 3, 4 } };

	std::vector< double > b{ 6, 7 };

	auto result = a * b;

	std::vector< double > correctResult{ 20, 46 };

	EXPECT_EQ( correctResult.size(), result.size() );

	for ( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6 );
	}

	/// Now test inverse multiplication
	result = b * a;

	std::vector< double > correctInverse{ 27, 40 };

	EXPECT_EQ( correctInverse.size(), result.size() );

	for ( auto i = 0u; i < correctInverse.size(); ++i ) {
		EXPECT_NEAR( correctInverse[ i ], result[ i ], 1e-6 );
	}
}