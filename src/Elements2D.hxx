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
		explicit ElementsLinear2D();

		FenestrationCommon::SquareMatrix< double > conductanceMatrix();
		// FenestrationCommon::SquareMatrix< double > & thermalCapacitanceMatrix();

		/// Creates lumped mass matrix that includes time derivative
		FenestrationCommon::Vector< double > getLumpedMass( const double DTime );
		FenestrationCommon::SquareMatrix< double > getMassMatrix( const double DTime );

		bool isLinear() const;

		void updateNodeValues( const std::vector< double > & values, const Property property );

		IElementLinear2D* findElement( const Node2D & t_Node1, const Node2D & t_Node2 );

	protected:
		/// FenestrationCommon::SquareMatrix< double > m_Conductance;
		/// FenestrationCommon::SquareMatrix< double > m_Capacitance;

		std::vector< std::unique_ptr< IElementLinear2D > > m_Elements;
		bool m_Linear;

	};

}