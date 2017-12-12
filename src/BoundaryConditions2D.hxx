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
    explicit BoundaryConditions2D( std::vector< std::unique_ptr< ILineLinear2D > > const & t_BCs );
		BoundaryConditions2D( const BoundaryConditions2D &other );
		BoundaryConditions2D& operator=( const BoundaryConditions2D& other );
    
    FenestrationCommon::CSquareMatrix* matrixA() const;
    std::vector< double > vectorR() const;
      
  private:
    std::unique_ptr< FenestrationCommon::CSquareMatrix > m_MatrixA;
    std::vector< double > m_vectorR;
  };
  
}