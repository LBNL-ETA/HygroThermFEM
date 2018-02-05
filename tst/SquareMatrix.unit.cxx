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
	SCOPED_TRACE( "Begin Test: Two elements example with transient." );
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
