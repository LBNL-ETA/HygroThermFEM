#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {
  
  class Elements2DLinear {
  public:
    explicit Elements2DLinear( const std::vector< ElementThermalLinear2D > & t_Elements );

    FenestrationCommon::SquareMatrix< double > & thermalConductivity();
    FenestrationCommon::SquareMatrix< double > & rhoCp();

  private:
    FenestrationCommon::SquareMatrix< double > m_Conductivity;
    FenestrationCommon::SquareMatrix< double > m_RhoCp;

  };
  
}