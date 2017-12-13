#pragma once

#include <memory>
#include "Element2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {
  
  class Elements2DLinear {
  public:
    explicit Elements2DLinear( const std::vector< ElementLinear2D > & t_Elements );

    FenestrationCommon::CSquareMatrix& thermalConductivity();
    FenestrationCommon::CSquareMatrix& rhoCp();

  private:
    FenestrationCommon::CSquareMatrix m_Conductivity;
    FenestrationCommon::CSquareMatrix m_RhoCp;

  };
  
}