#pragma once

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"

namespace MoisThermFEM {

	/// Class to keep elements and boundary conditions together. One domain per Thermal, Mass and
	/// Pressure governing equations. It makes system of equations necessary to be solved
	class Domain {
	public:
		Domain( const ElementsLinear2D & m_Elements, const BoundaryConditions2D & m_BCs );		

		std::vector< double > steadyState();

		std::vector< double > transient( std::vector< double > & currentTimestepValues, const double t_DTime );

	protected:

		FenestrationCommon::SquareMatrix< double > steadyStateLeftHandSide();
		std::vector< double > steadyStateRightHandSide();

		FenestrationCommon::SquareMatrix< double > transientLeftHandSide( const double t_DTime );
		std::vector< double > transientRightHandSide( std::vector< double > & t_PreviousSolution, const double t_DTime );

		ElementsLinear2D m_Elements;
		BoundaryConditions2D m_BCs;
	};

}