#include <cmath>
#include <algorithm>

#include "Domains.hxx"
#include "NodePool.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"

using FenestrationCommon::CLinearSolver;

namespace MoisThermFEM {

	Domain::Domain( const ElementsLinear2D & t_Elements, const BoundaryConditions2D & t_BCs )
			: m_Elements( t_Elements ), m_BCs( t_BCs ), m_Linear( t_BCs.isLinear() ) {}

	FenestrationCommon::SquareMatrix< double > Domain::steadyStateLeftHandSide() {
		auto condMat = m_Elements.conductanceMatrix();
		auto H = m_BCs.HMatrix();
		condMat = condMat +  H;

		return condMat;
	}

	std::vector< double > Domain::steadyStateRightHandSide() {
		return m_BCs.RVector();
	}

	FenestrationCommon::SquareMatrix< double > Domain::transientM_K_H_Matrix( const double t_DTime ) {
		std::vector< double > M = m_Elements.getLumpedMass( t_DTime );
		auto conductanceMatrix = m_Elements.conductanceMatrix();
		conductanceMatrix = conductanceMatrix.addDiagonal( M );
		conductanceMatrix = conductanceMatrix + m_BCs.HMatrix();

		return conductanceMatrix;
	}

	std::vector< double > Domain::transientMT_R_Vector( std::vector< double > & t_PreviousSolution,
																											const double t_DTime ) {
		std::vector< double > M = m_Elements.getLumpedMass( t_DTime );
		std::vector< double > B( t_PreviousSolution.size() );

		auto size = NodePool::Instance().maxIndex();
		auto Rs = m_BCs.RVector();

		for ( auto j = 0u; j < size; ++j ) {
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

		auto B = transientMT_R_Vector( currentTimestepValues, t_DTime );
		auto A = transientM_K_H_Matrix( t_DTime );

		CLinearSolver aSolver;

		std::vector< double > solution;

		if( m_Linear ) {
			solution = aSolver.solveSystem( A, B );
		} else {
			solution = currentTimestepValues;

			auto error = std::numeric_limits< double >::max();

			size_t numOfIterations = 0;

			while( error > ConvergenceError ) {
				auto temp = A * solution;
				/// temp = B - temp
				std::transform( B.begin(), B.end(), temp.begin(),
												temp.begin(), std::minus< double >() );

				auto dU = aSolver.solveSystem( A, temp );

				error = norm( dU );

				std::transform( dU.begin(), dU.end(),
												solution.begin(),
												solution.begin(), std::plus< double >() );
				++numOfIterations;
				if( numOfIterations > MaxIterations ) {
					throw std::runtime_error("Solution failed to converge.");
				}
			}

		}

		return solution;
	}

	double Domain::norm( const std::vector< double > & t_vector ) {
		double result{ 0 };
		std::for_each( t_vector.begin(), t_vector.end(), [ & ]( double n ) {
			result += n * n;
		} );

		result = std::pow( result, 0.5 );

		return result;
	}
}