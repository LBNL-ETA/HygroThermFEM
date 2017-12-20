#include <memory>
#include <gtest/gtest.h>

#include "Conrad2D.hxx"

using namespace FenestrationCommon;

class TestLinearSolver1 : public testing::Test {

private:
	CLinearSolver m_Solver;

protected:
	void
	SetUp() override
	{
		
	}

public:
	CLinearSolver&
	GetSolver()	{
		return m_Solver;
	}

};

TEST_F( TestLinearSolver1, Test1 )
{
	SCOPED_TRACE( "Begin Test: Test Linear Solver (1) - Solving simple matrix." );

	const auto size = 3;
	SquareMatrix< double > aMatrix( size );

	aMatrix[ 0 ] = { 2, 1, 3 };
	aMatrix[ 1 ] = { 2, 6, 8 };
	aMatrix[ 2 ] = { 6, 8, 18 };

	std::vector< double > aVector = { 1, 3, 5 };

	auto aSolver = GetSolver();

	auto aSolution = aSolver.solveSystem( aMatrix, aVector );

	EXPECT_NEAR( 3.0/10.0, aSolution[ 0 ], 1e-6 );
	EXPECT_NEAR( 2.0/5.0, aSolution[ 1 ], 1e-6 );
	EXPECT_NEAR( 0.0, aSolution[ 2 ], 1e-6 );
}

TEST_F( TestLinearSolver1, Test2 )
{
	SCOPED_TRACE( "Begin Test: Test Linear Solver (2) - Solving simple matrix." );

	const auto size = 4;
	SquareMatrix< double > aMatrix( size );

	// This set-up is from fortran tarcog example
	aMatrix[ 0 ] = { 32817.2867004354, 1, 0, -32808.3972386696 };
	aMatrix[ 1 ] = { 1.28054053432588, -1, 0, 0 };
	aMatrix[ 2 ] = { 0, 0, -1, 1.26433319889839 };
	aMatrix[ 3 ] = { 32808.3972386696, 0, -1, -32810.4664383299 };

	std::vector< double > aVector = { 3163.241853, -73.479324, -67.913411, -1070.271453 };

	auto aSolver = GetSolver();

	auto aSolution = aSolver.solveSystem( aMatrix, aVector );

	EXPECT_NEAR( 303.040746, aSolution[ 0 ], 1e-6 );
	EXPECT_NEAR( 461.535283, aSolution[ 1 ], 1e-6 );
	EXPECT_NEAR( 451.057585, aSolution[ 2 ], 1e-6 );
	EXPECT_NEAR( 303.040507, aSolution[ 3 ], 1e-6 );
}
