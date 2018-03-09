#include <memory>
#include <gtest/gtest.h>

#include "MoisThermFEM2D.hxx"

using FenestrationCommon::Vector;

class VectorTest : public testing::Test {

protected:
	void
	SetUp() override {
	}

	void
	TearDown() override {

	}

};

TEST_F( VectorTest, TestAddition ) {
	SCOPED_TRACE( "Begin Test: Vector addition." );
	const Vector< double > a{ 1, 2, 3 };
	const Vector< double > b{ 4, 5, 6 };

	auto result = a + b;

	Vector< double > correctResult{ 5, 7, 9 };

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6  );
	}


	/// Now test commutation
	result = b + a;

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6 );
	}

}

TEST_F( VectorTest, TestSubraction ) {
	SCOPED_TRACE( "Begin Test: Vector subtraction." );

	const Vector< double > a{ 1, 2, 3 };
	const Vector< double > b{ 4, 5, 6 };

	auto result = a - b;

	Vector< double > correctResult{ -3, -3, -3 };

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6  );
	}


	/// Now test commutation
	result = b - a;

	Vector< double > correctInverse{ 3, 3, 3 };

	EXPECT_EQ( result.size(), correctInverse.size() );

	for( auto i = 0u; i < correctInverse.size(); ++i ) {
		EXPECT_NEAR( correctInverse[ i ], result[ i ], 1e-6 );
	}
	
}

TEST_F( VectorTest, TestMultiplication ) {
	SCOPED_TRACE( "Begin Test: Vector." );
	const Vector< double > a{ 1, 2, 3 };
	const Vector< double > b{ 4, 5, 6 };

	auto result = a * b;

	Vector< double > correctResult{ 4, 10, 18 };

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6  );
	}


	/// Now test commutation
	result = b * a;

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6 );
	}
}

TEST_F( VectorTest, TestMultiplicationWithSingleValue ) {
	SCOPED_TRACE( "Begin Test: Vector multiplication with single value." );

	const Vector< double > a{ 1, 2, 3 };
	const double b{ 4 };

	auto result = a * b;

	Vector< double > correctResult{ 4, 8, 12 };

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6  );
	}

	/// Now test commutation
	result = b * a;

	EXPECT_EQ( result.size(), correctResult.size() );

	for( auto i = 0u; i < correctResult.size(); ++i ) {
		EXPECT_NEAR( correctResult[ i ], result[ i ], 1e-6 );
	}
	
}