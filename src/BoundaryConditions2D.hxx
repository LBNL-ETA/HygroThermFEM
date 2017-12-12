#pragma once

#include <vector>
#include <memory>

namespace FenestrationCommon {
  class CSquareMatrix;
}

namespace MoisThermFEM {
  
  class ILineLinear2D;
  
  class BoundaryConditions2D {
  public:
    explicit BoundaryConditions2D( std::vector< std::shared_ptr< ILineLinear2D > > const & t_BCs );
    
    std::shared_ptr< FenestrationCommon::CSquareMatrix > matrixA() const;
    std::vector< double > vectorR() const;
      
  private:
    std::shared_ptr< FenestrationCommon::CSquareMatrix > m_MatrixA;
    std::vector< double > m_vectorR;
  };
  
}