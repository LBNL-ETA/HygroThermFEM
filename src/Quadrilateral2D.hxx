#ifndef QUADRILATERAL2D_H_9ea233f0_7531_4750_a6fc_eab6807003bf
#define QUADRILATERAL2D_H_9ea233f0_7531_4750_a6fc_eab6807003bf

#include <vector>
#include <memory>

#include "Node2D.hxx"

namespace Conrad {

  ////////////////////////////////////////////////////////////////////////////
  //   QuadrilateralLinearGlobal2D
  ////////////////////////////////////////////////////////////////////////////
  class QuadrilateralLinearGlobal2D {
  public:
    QuadrilateralLinearGlobal2D( 
      Node2D const & t_Node1,
      Node2D const & t_Node2,
      Node2D const & t_Node3,
      Node2D const & t_Node4 );

    QuadrilateralLinearGlobal2D( const QuadrilateralLinearGlobal2D& t_Element );

    // Coordinates of integration (Gauss) point in global coordinate system
    double xg( size_t const IntPointIndex ) const;
    double yg( size_t const IntPointIndex ) const;

    // Shape function derivatives in global coordinate system
    std::vector< double > DPsiDx( size_t const IntPointIndex ) const;
    std::vector< double > DPsiDy( size_t const IntPointIndex ) const;

    double det( size_t const IntPointIndex ) const;

    // Element nodes
    std::vector< Node2D > getNodes() const;

  private:
    QuadrilateralNodes2D m_Nodes;

    ////////////////////////////////////////////////////////////////////////////
    //   GaussPoint2DGlobal
    ////////////////////////////////////////////////////////////////////////////

    // Handles single integration point in global coordinate space (cartesian 2D)
    class GaussPoint2DGlobal {
    public:
      // Nodes represent global coordinates and Index is the index of integration point
      GaussPoint2DGlobal( 
        Node2D const & t_Node1,
        Node2D const & t_Node2,
        Node2D const & t_Node3,
        Node2D const & t_Node4,
        size_t const Index );

      size_t getIndex() const;

      double xg() const;
      double yg() const;

      std::vector< double > getDPsiDx() const;
      std::vector< double > getDPsiDy() const;

      double det() const;

    private:
      size_t m_Index; // Index of gauss point withing global element
      double m_Xg;
      double m_Yg;
      double m_JacobiDet;
      std::vector< double > m_DPsiDx;
      std::vector< double > m_DPsiDy;

    };

    std::vector< GaussPoint2DGlobal > m_GaussPoints;

  };
  
}

#endif