#include <assert.h>

#include "Node2D.hxx"

namespace Conrad {

  ////////////////////////////////////////////////////////////////////////////
  //   LocalPoint1D
  ////////////////////////////////////////////////////////////////////////////

  LocalPoint1D::LocalPoint1D( double const t_ksi ) : ksi( t_ksi ) {

  }

  LocalPoint1D::LocalPoint1D( LocalPoint1D const & t_LocalPoint ) {
    ksi = t_LocalPoint.ksi;
  }

  ////////////////////////////////////////////////////////////////////////////
  //   LocalPoint2D
  ////////////////////////////////////////////////////////////////////////////

  LocalPoint2D::LocalPoint2D( double const t_ksi, double const t_eta ) :
    ksi( t_ksi ), eta( t_eta ) {

  };
  
  LocalPoint2D::LocalPoint2D( LocalPoint2D const & t_LocalPoint ) {
    ksi = t_LocalPoint.ksi;
    eta = t_LocalPoint.eta;
  };

  ////////////////////////////////////////////////////////////////////////////
  //   Node2D
  ////////////////////////////////////////////////////////////////////////////
  
  Node2D::Node2D( 
    size_t const t_NodeNumber, 
    double const t_x, 
    double const t_y, 
    double const t_temperature ) :
    nodeNumber( t_NodeNumber ), x( t_x ), y( t_y ), temperature( t_temperature ) {
    
  };

  Node2D::Node2D( Node2D const & t_Node ) {
    nodeNumber = t_Node.nodeNumber;
    x = t_Node.x;
    y = t_Node.y;
    temperature = t_Node.temperature;
  };

  ////////////////////////////////////////////////////////////////////////////
  //   INodesStorage
  ////////////////////////////////////////////////////////////////////////////

  INodesStorage::INodesStorage() {

  }

  Node2D INodesStorage::getNode( size_t const Index ) const {
    assert( Index < m_Nodes.size() );
    return m_Nodes[ Index ];
  }

  std::vector< size_t > INodesStorage::getNodeIndexes() const {
    std::vector< size_t > indexes;
    for( const Node2D& aNode : m_Nodes ) {
      indexes.push_back( aNode.nodeNumber );
    }
    return indexes;
  }

  ////////////////////////////////////////////////////////////////////////////
  //   LineNodes2D
  ////////////////////////////////////////////////////////////////////////////

  LineNodes2D::LineNodes2D( 
    Node2D const & t_Node1, 
    Node2D const & t_Node2 ) : INodesStorage() {
    m_Nodes.push_back( t_Node1 );
    m_Nodes.push_back( t_Node2 );
  }

  ////////////////////////////////////////////////////////////////////////////
  //   QuadrilateralNodes2D
  ////////////////////////////////////////////////////////////////////////////

  QuadrilateralNodes2D::QuadrilateralNodes2D( 
    Node2D const & t_Node1, 
    Node2D const & t_Node2, 
    Node2D const & t_Node3, 
    Node2D const & t_Node4 ) {
    m_Nodes.push_back( t_Node1 );
    m_Nodes.push_back( t_Node2 );
    m_Nodes.push_back( t_Node3 );
    m_Nodes.push_back( t_Node4 );
  }

}