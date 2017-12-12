#pragma once

#include <memory>
#include "Element2D.hxx"

namespace FenestrationCommon {
  class CSquareMatrix;
}

namespace MoisThermFEM {
  
  class Elements2DLinear {
  public:
    explicit Elements2DLinear( std::vector< ElementLinear2D > const & t_Elements );

    std::shared_ptr< FenestrationCommon::CSquareMatrix > thermalConductivity() const;
    std::shared_ptr< FenestrationCommon::CSquareMatrix > rhoCp() const;

  private:
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_Conductivity;
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_RhoCp;

  };
  
}