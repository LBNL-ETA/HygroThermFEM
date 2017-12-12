#pragma once

#include <vector>

#include "Node2D.hxx"

namespace MoisThermFEM {

  struct Node2D;
  class LineNodes2D;
  class LineLinearLocal1D;
  
  // Interface class for all boundary conditions
  class ILineLinear2D {
  public:
    ILineLinear2D( 
      Node2D const & t_Node1, 
      Node2D const & t_Node2 );

    std::vector< std::size_t > getNodeIndexes() const;
    std::vector< double > rightHandSideVector() const;
    std::vector< std::vector< double > > matrixA() const;

  protected:
    static std::size_t numOfIntegrationPoints();
    static double Psi( std::size_t const IntegrationPointIndex, std::size_t const Index );

    LineNodes2D m_Nodes;
    double m_Determinant;
    std::vector< double > m_Rvector; // right hand side vector
    std::vector< std::vector< double > > m_matrixA;

  };
  
}