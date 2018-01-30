#include "Domains.hxx"
#include "NodePool.hxx"
#include "LinearSolver.hxx"

using FenestrationCommon::CLinearSolver;

namespace MoisThermFEM {

	Domain::Domain( const ElementsLinear2D & m_Elements, const BoundaryConditions2D & m_BCs )
			: m_Elements( m_Elements ), m_BCs( m_BCs ) {}

	FenestrationCommon::SquareMatrix< double > Domain::steadyStateLeftHandSide() {
		auto condMat = m_Elements.thermalConductanceMatrix();
		auto H = m_BCs.matrixA();
		condMat = condMat.add( H );

		return condMat;
	}

	std::vector< double > Domain::steadyStateRightHandSide() {
		return m_BCs.vectorR();
	}

	FenestrationCommon::SquareMatrix< double > Domain::transientLeftHandSide( const double t_DTime ) {
		std::vector< double > M = m_Elements.getLumpedMass( t_DTime );
		auto conductanceMatrix = m_Elements.thermalConductanceMatrix();
		conductanceMatrix = conductanceMatrix.addDiagonal( M );
		conductanceMatrix = conductanceMatrix.add( m_BCs.matrixA() );

		return conductanceMatrix;
	}

	std::vector< double > Domain::transientRightHandSide( std::vector< double > & t_PreviousSolution,
																												const double t_DTime ) {
		std::vector< double > M = m_Elements.getLumpedMass( t_DTime );
		std::vector< double > B( t_PreviousSolution.size() );

		auto size = NodePool::Instance().maxIndex();
		auto Rs = m_BCs.vectorR();

		for ( unsigned j = 0; j < size; ++j ) {
			B[ j ] = t_PreviousSolution[ j ] * M[ j ] + Rs[ j ];
		}

		return B;
	}

	std::vector< double > Domain::steadyState() {
		auto B = steadyStateRightHandSide();

		CLinearSolver aSolver;

		return aSolver.solveSystem( steadyStateLeftHandSide(), B );
	}

	std::vector< double >
	Domain::transient( std::vector< double > & currentTimestepValues, const double t_DTime ) {

		auto B = transientRightHandSide( currentTimestepValues, t_DTime );
		auto A = transientLeftHandSide( t_DTime );

		CLinearSolver aSolver;

		return aSolver.solveSystem( A, B );
	}
}