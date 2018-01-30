#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {
  
	// Container class to hold onto thermal conductanceMatrix and heat capacity for
	// all nodes in the domain.
  class ElementsThermalLinear2D {
  public:
    explicit ElementsThermalLinear2D( const std::vector< ElementThermalLinear2D > & t_Elements );

    FenestrationCommon::SquareMatrix< double > & conductanceMatrix();
    FenestrationCommon::SquareMatrix< double > & capacitanceMatrix();

  private:
    FenestrationCommon::SquareMatrix< double > m_Conductance;
    FenestrationCommon::SquareMatrix< double > m_Capacitance;

  };
  
}