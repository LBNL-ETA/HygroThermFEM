#include <math.h>
#include <assert.h>

#include "Line2D.hxx"
#include "Node2D.hxx"
#include "IntegrationPoints.hxx"
#include "LineLocal1D.hxx"

namespace Conrad {

  ILineLinear2D::ILineLinear2D( 
    Node2D const & t_Node1,
    Node2D const & t_Node2 ) : m_Nodes( t_Node1, t_Node2 ) {

    m_Determinant = 0.5 * sqrt( pow( t_Node1.x - t_Node2.x, 2 ) + pow( t_Node1.y - t_Node2.y, 2 ) );

    // This will just reserve correct room for vectors and matrices
    m_Rvector.resize( 2 );
    m_matrixA.resize( 2 );
    for( auto i = 0u; i < m_matrixA.size(); ++i ) {
      m_matrixA[ i ].resize( 2 );
    }
  }

  std::vector< size_t > ILineLinear2D::getNodeIndexes() const {
    return m_Nodes.getNodeIndexes();
  }

  std::vector< double > ILineLinear2D::rVector() const {
    return m_Rvector;
  }

  std::vector< std::vector< double > > ILineLinear2D::matrixA() const {
    return m_matrixA;
  }

  size_t ILineLinear2D::numOfIntegrationPoints() {
    return IntegrationPoints2D::Instance().count1D();
  }

  double ILineLinear2D::Psi( size_t const IntegrationPointIndex, size_t const Index ) {
    return LineLinearLocal1D::Instance().Psi( IntegrationPointIndex, Index );
  }

}