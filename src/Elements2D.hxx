#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {
  
	// Container class to hold onto thermal conductivity and heat capacity for 
	// all nodes in the domain.
  class ElementsThermalLinear2D {
  public:
    explicit ElementsThermalLinear2D( const std::vector< ElementThermalLinear2D > & t_Elements );

    FenestrationCommon::SquareMatrix< double > & thermalConductivity();
    FenestrationCommon::SquareMatrix< double > & rhoCp();

  private:
    FenestrationCommon::SquareMatrix< double > m_Conductivity;
    FenestrationCommon::SquareMatrix< double > m_RhoCp;

  };
  
}