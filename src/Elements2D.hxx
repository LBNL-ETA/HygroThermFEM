#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {
  
	// Container class to hold all elements connected into global matrix
  class ElementsLinear2D {
  public:
    explicit ElementsLinear2D( const std::vector< ElementLinear2D > & t_Elements );

    FenestrationCommon::SquareMatrix< double > & thermalConductanceMatrix();
    FenestrationCommon::SquareMatrix< double > & thermalCapacitanceMatrix();

		/// Creates lumped mass matrix that includes time derivative
		std::vector< double > getLumpedMass( const double DTime );

  private:
    FenestrationCommon::SquareMatrix< double > m_Conductance;
    FenestrationCommon::SquareMatrix< double > m_Capacitance;

  };
  
}