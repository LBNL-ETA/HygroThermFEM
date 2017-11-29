#include "Elements2D.hxx"
#include "NodePool.hxx"
#include "SquareMatrix.hxx"

using namespace FenestrationCommon;

namespace Conrad {

  Elements2DLinear::Elements2DLinear( std::vector< ElementLinear2D > const & t_Elements ) {

    // Size is kept in nood pool
    auto maxIndex = NodePool::Instance().maxIndex();

    m_Conductivity = std::make_shared< CSquareMatrix >( maxIndex );
    m_RhoCp = std::make_shared< CSquareMatrix >( maxIndex );

    // now copy matrices into global matrix
    for( auto const & aElement : t_Elements ) {
      auto indexes = aElement.nodeIndexes();
      auto conductivity = aElement.conductivity();
      auto rhoCp = aElement.RhoCp();
      for( size_t i = 0; i < 4; ++i ) {
        for( size_t j = 0; j < 4; ++j ) {
          ( *m_Conductivity )[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += conductivity[ i ][ j ];
          ( *m_RhoCp )[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += rhoCp[ i ][ j ];
        }
      }
    }

  }

  std::shared_ptr< CSquareMatrix > Elements2DLinear::conductivity() const {
    return m_Conductivity;
  }

  std::shared_ptr< CSquareMatrix > Elements2DLinear::rhoCp() const {
    return m_RhoCp;
  }

}