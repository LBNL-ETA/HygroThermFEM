#pragma once

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"

namespace MoisThermFEM {

	/// Class to keep elements and boundary conditions together. One domain per Thermal, Mass and
	/// Pressure governing equations.
	class Domain {
	public:
		Domain( const ElementsLinear2D & m_Elements, const BoundaryConditions2D & m_BCs );

		/// Calculates steady state solution
		std::vector< double > steadyState();

		/// Calculates next timestep value
		std::vector< double >
		transient( std::vector< double > & currentTimestepValues, const double t_DTime );

	protected:

		FenestrationCommon::SquareMatrix< double > steadyStateLeftHandSide();
		std::vector< double > steadyStateRightHandSide();

		/// In matrix equations some structures are showing up in both (linear and nonlinear) cases
		/// and those matrix operations are separated into functions.
		/// This function retrieves M+K+H matrix
		FenestrationCommon::SquareMatrix< double > transientM_K_H_Matrix( const double t_DTime );

		/// This function retrieves M*U+R vector (where U is state variable)
		std::vector< double >
		transientMT_R_Vector( std::vector< double > & t_PreviousSolution, const double t_DTime );

		/// Returns norm of vector. Necessary to estimate current error.
		static double norm( const std::vector< double > & t_vector );

		ElementsLinear2D m_Elements;
		BoundaryConditions2D m_BCs;
		bool m_Linear;
	};

}