#pragma once

#include <vector>

namespace Conrad {

  ////////////////////////////////////////////////////////////////////////////
  //   LocalPoint1D
  ////////////////////////////////////////////////////////////////////////////

  struct LocalPoint1D {
    explicit LocalPoint1D( double const t_ksi );
    LocalPoint1D( LocalPoint1D const & t_LocalPoint );
    double ksi;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   LocalPoint2D
  ////////////////////////////////////////////////////////////////////////////

  struct LocalPoint2D {
    LocalPoint2D( double const t_ksi, double const t_eta );
    LocalPoint2D( LocalPoint2D const & t_LocalPoint );
    double ksi;
    double eta;

  };


  ////////////////////////////////////////////////////////////////////////////
  //   Node2D
  ////////////////////////////////////////////////////////////////////////////

  // Defines nodal point in two dimensional cartesian space
  struct Node2D {
    Node2D( std::size_t const t_NodeNumber, double const t_x, double const t_y, double const y_temperature = 0 );
    Node2D( Node2D const & t_Node );

    std::size_t nodeNumber;
    double x;
    double y;
    double temperature;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   INodesStorage
  ////////////////////////////////////////////////////////////////////////////
  class INodesStorage {
  public:
    INodesStorage();

    Node2D getNode( size_t const Index ) const;
    std::vector< size_t > getNodeIndexes() const;

  protected:
    std::vector< Node2D > m_Nodes;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   LineNodes2D
  ////////////////////////////////////////////////////////////////////////////
  class LineNodes2D : public INodesStorage {
  public:
    LineNodes2D( 
      Node2D const & t_Node1, 
      Node2D const & t_Node2 );

  };

  ////////////////////////////////////////////////////////////////////////////
  //   QuadrilateralNodes2D
  ////////////////////////////////////////////////////////////////////////////
  class QuadrilateralNodes2D : public INodesStorage {
  public:
    QuadrilateralNodes2D(
      Node2D const & t_Node1, 
      Node2D const & t_Node2, 
      Node2D const & t_Node3, 
      Node2D const & t_Node4 );
  };
  
}