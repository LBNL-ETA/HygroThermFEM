#pragma once

#include <memory>
#include <functional>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {

	/// Container class to hold all elements connected into global matrix. This is only for elements
	/// and not for boundary conditions
	class ElementsLinear2D {
	public:
		explicit ElementsLinear2D(
				const std::vector< std::reference_wrapper< const IElementLinear2D > > & t_Elements );

		FenestrationCommon::SquareMatrix< double > & conductanceMatrix();
		// FenestrationCommon::SquareMatrix< double > & thermalCapacitanceMatrix();

		/// Creates lumped mass matrix that includes time derivative
		FenestrationCommon::Vector< double > getLumpedMass( const double DTime );

	private:
		FenestrationCommon::SquareMatrix< double > m_Conductance;
		FenestrationCommon::SquareMatrix< double > m_Capacitance;

	};

}