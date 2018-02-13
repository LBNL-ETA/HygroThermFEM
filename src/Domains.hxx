#pragma once

#include "Elements2D.hxx"
#include "BoundaryCondition2DFactory.hxx"
#include "Element2DFactory.hxx"

namespace MoisThermFEM {

	/// Class to keep elements and boundary conditions together. One domain per Thermal, Mass and
	/// Pressure governing equations.
	class Domain {
	public:
		Domain() = default;

		/// Calculates steady state solution
		std::vector< double > steadyState();

		/// Calculates next timestep value from current values
		std::vector< double >
		transient( std::vector< double > & currentStateValues, const double t_DTime );

		BoundaryCondition2DFactory & boundariesCreator();
		Element2DFactory & elementsCreator();

	protected:

		FenestrationCommon::SquareMatrix< double > steadyStateLeftHandSide();
		FenestrationCommon::Vector< double > steadyStateRightHandSide();

		/// In matrix equations some structures are showing up in both (linear and nonlinear) cases
		/// and those matrix operations are separated into functions.
		/// This function retrieves M+K+H matrix
		FenestrationCommon::SquareMatrix< double > transientM_K_H_Matrix( const double t_DTime );
		/// FenestrationCommon::SquareMatrix< double > transientDH_Matrix();

		/// This function retrieves M*U+R vector (where U is state variable)
		FenestrationCommon::Vector< double >
		transientMT_R_Vector( std::vector< double > & t_PreviousSolution, const double t_DTime );

		/// Returns norm of vector. Necessary to estimate current error.
		static double norm( const std::vector< double > & t_vector );

		bool isLinear() const;

		Element2DFactory m_Elements;
		BoundaryCondition2DFactory m_BCs;
	};

}