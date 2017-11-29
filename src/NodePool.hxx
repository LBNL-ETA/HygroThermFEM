#pragma once

#include <vector>

namespace Conrad {
  
  struct Node2D;
  
  class NodePool {
  public:
    static NodePool& Instance();
    
    Node2D& createNode(
      size_t const t_NodeNumber,
      double const t_x,
      double const t_y,
      double const t_temperature = 0 );

    Node2D& getNode( size_t const Index );
    size_t maxIndex() const;
    std::vector< double > nodeTemperatures() const;
    void clear();
    
  private:
    NodePool();
    ~NodePool();
    
    std::vector< Node2D > m_Nodes;
  };
  
}