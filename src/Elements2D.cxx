#include "Elements2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

  Elements2DLinear::Elements2DLinear( std::vector< ElementLinear2D > const & t_Elements ) : 
		m_Conductivity( NodePool::Instance().maxIndex() ), m_RhoCp( NodePool::Instance().maxIndex() ) {

    // now copy matrices into global matrix
    for( auto const & aElement : t_Elements ) {
      auto indexes = aElement.nodeIndexes();
      auto conductivity = aElement.conductivity();
      auto rhoCp = aElement.RhoCp();
      for( size_t i = 0; i < 4; ++i ) {
        for( size_t j = 0; j < 4; ++j ) {
          m_Conductivity[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += conductivity[ i ][ j ];
          m_RhoCp[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += rhoCp[ i ][ j ];
        }
      }
    }

  }

  CSquareMatrix& Elements2DLinear::thermalConductivity() {
    return m_Conductivity;
  }

  CSquareMatrix& Elements2DLinear::rhoCp() {
    return m_RhoCp;
  }

}