#include <cmath>
#include <algorithm>

#include "Domains.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"

using FenestrationCommon::CLinearSolver;

namespace MoisThermFEM {

	Domain::Domain( const Property property ) : m_Property( property ) {

	}

	FenestrationCommon::SquareMatrix< double > Domain::steadyStateLeftHandSide() {
		auto condMat = m_Elements.conductanceMatrix();
		auto H = m_BCs.HMatrix();
		condMat = condMat + H;

		return condMat;
	}

	FenestrationCommon::Vector< double > Domain::steadyStateRightHandSide() {
		return m_BCs.RVector();
	}

	FenestrationCommon::SquareMatrix< double > Domain::transientM_K_H_Matrix( const double t_DTime ) {
		auto M = m_Elements.getLumpedMass( t_DTime );
		// auto M = m_Elements.getMassMatrix( t_DTime );
		auto M_K_H = m_Elements.conductanceMatrix();
		M_K_H = M_K_H.addDiagonal( M );
		// M_K_H += M;
		M_K_H += m_BCs.HMatrix();

		return M_K_H;
	}

	FenestrationCommon::Vector< double >
	Domain::transientMT_R_Vector( std::vector< double > & t_PreviousSolution,
																const double t_DTime ) {
		FenestrationCommon::Vector< double > M{ m_Elements.getLumpedMass( t_DTime ) };
		auto R = m_BCs.RVector();

		auto B = t_PreviousSolution * M + R;

		return B;
	}

	std::vector< double > Domain::steadyState() {
		auto B = steadyStateRightHandSide();
		CLinearSolver aSolver;
		return aSolver.solveSystem( steadyStateLeftHandSide(), B );
	}

	std::vector< double >
	Domain::transient( std::vector< double > & currentStateValues, const double t_DTime ) {

		auto A = transientM_K_H_Matrix( t_DTime );
		auto B = transientMT_R_Vector( currentStateValues, t_DTime );

		CLinearSolver aSolver;

		std::vector< double > solution;

		if( isLinear() ) {
			solution = aSolver.solveSystem( A, B );
		} else {
			solution = currentStateValues;

			auto error = std::numeric_limits< double >::max();

			size_t numOfIterations = 0;

			while ( error > ConvergenceError ) {
				auto temp = A * solution;
				temp = B - temp;

				/// Seems that DH can be avoided. Same solution is achieved faster without it. Topaz does
				/// have this implementation. Will keep it commented in case we want to test it in future
				/// when new kind of boundary conditions are introduced (Simon)
				/// auto DH = transientDH_Matrix( );
				/// DH = A + DH;

				/// auto dU = aSolver.solveSystem( DH, temp );

				auto dU = aSolver.solveSystem( A, temp );

				error = norm( dU );

				std::transform( dU.begin(), dU.end(),
												solution.begin(),
												solution.begin(), std::plus< double >() );

				++numOfIterations;

				m_BCs.updateNodeValues( solution, m_Property );

				A = transientM_K_H_Matrix( t_DTime );
				B = transientMT_R_Vector( currentStateValues, t_DTime );

				if( numOfIterations > MaxIterations ) {
					throw std::runtime_error( "Solution failed to converge." );
				}
			}

		}

		return solution;
	}

	BoundaryCondition2DFactory & Domain::boundariesCreator() {
		return m_BCs;
	}

	Element2DFactory & Domain::elementsCreator() {
		return m_Elements;
	}

	bool Domain::isLinear() const {
		return m_BCs.isLinear() && m_Elements.isLinear();
	}

	void Domain::updateNodeValues( const std::vector< double > & values, const Property property ) {
		m_BCs.updateNodeValues( values, property );
		m_Elements.updateNodeValues( values, property );
	}

	/// FenestrationCommon::SquareMatrix< double > Domain::transientDH_Matrix() {
	/// 	return m_BCs.DHMatrix();;
	/// }
}