#include "BoundaryConditions2D.hxx"
#include "Line2D.hxx"
#include "NodePool.hxx"
#include "SquareMatrix.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {
  
  BoundaryConditions2D::BoundaryConditions2D( std::vector< std::shared_ptr< ILineLinear2D > > const & t_BCs ) {
    
    // Size is kept in node pool
    auto maxIndex = NodePool::Instance().maxIndex();

    m_vectorR.resize( maxIndex );
    m_MatrixA = std::make_shared< CSquareMatrix >( maxIndex );

    // Create full size matrices
    for( const auto& aBc : t_BCs ) {
      auto indexes = aBc->getNodeIndexes();
      auto matA = aBc->matrixA();
      auto vecR = aBc->rightHandSideVector();
      for( size_t i = 0; i < 2; ++i ) {
        for( size_t j = 0; j < 2; ++j ) {
          ( *m_MatrixA )[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += matA[ i ][ j ];
        }
        m_vectorR[ indexes[ i ] - 1 ] += vecR[ i ];
      }
    }
  }

  std::shared_ptr< CSquareMatrix > BoundaryConditions2D::matrixA() const {
    return m_MatrixA;
  }

  std::vector< double > BoundaryConditions2D::vectorR() const {
    return m_vectorR;
  }

}