#include <assert.h>

#include "Node2D.hxx"

namespace MoisThermFEM {

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
  }
  
  LocalPoint2D::LocalPoint2D( LocalPoint2D const & t_LocalPoint ) {
    ksi = t_LocalPoint.ksi;
    eta = t_LocalPoint.eta;
  }

  ////////////////////////////////////////////////////////////////////////////
  //   Node2D
  ////////////////////////////////////////////////////////////////////////////
  
  Node2D::Node2D( 
    const std::size_t t_NodeNumber,
    const double t_x, 
    const double t_y,
		const double t_temperature ) : 
		nodeNumber( t_NodeNumber ), x( t_x ), y( t_y ), temperature( t_temperature ) {
    
  }

  Node2D::Node2D( const Node2D & t_Node ) {
    nodeNumber = t_Node.nodeNumber;
    x = t_Node.x;
    y = t_Node.y;
		temperature = t_Node.temperature;
  }


  ////////////////////////////////////////////////////////////////////////////
  //   INodesStorage
  ////////////////////////////////////////////////////////////////////////////

  INodesStorage::INodesStorage() {

  }

	INodesStorage::INodesStorage( std::initializer_list< Node2D > t_Nodes ) : m_Nodes( t_Nodes ) {
	}

	Node2D INodesStorage::getNode( const std::size_t Index ) const {
    assert( Index < m_Nodes.size() );
    return m_Nodes[ Index ];
  }

  std::vector< std::size_t > INodesStorage::getNodeIndexes() const {
    std::vector< std::size_t > indexes;
    for( const auto & aNode : m_Nodes ) {
      indexes.push_back( aNode.nodeNumber );
    }
    return indexes;
  }

  ////////////////////////////////////////////////////////////////////////////
  //   LineNodes2D
  ////////////////////////////////////////////////////////////////////////////

  LineNodes2D::LineNodes2D( 
    const Node2D & t_Node1, 
    const Node2D & t_Node2 ) : INodesStorage( { t_Node1, t_Node2 } ) {
    
  }

  ////////////////////////////////////////////////////////////////////////////
  //   QuadrilateralNodes2D
  ////////////////////////////////////////////////////////////////////////////

  QuadrilateralNodes2D::QuadrilateralNodes2D( 
    const Node2D & t_Node1, 
    const Node2D & t_Node2, 
    const Node2D & t_Node3, 
		const Node2D & t_Node4 ) : INodesStorage( { t_Node1, t_Node2, t_Node3, t_Node4 } ) {
    
  }

}