#include <cassert>
#include <algorithm>

#include "NodePool.hxx"

#ifdef STL_MULTITHREADING
#    include <execution>
#endif

namespace HygroThermFEM
{
    NodePool & NodePool::Instance()
    {
        static NodePool m_Instance;
        return m_Instance;
    }

    Node2D & NodePool::createNode(const std::size_t t_NodeNumber,
                                  const double t_x,
                                  const double t_y,
                                  const State & t_Prop)
    {
        m_Nodes.emplace_back(t_NodeNumber, t_x, t_y, t_Prop);
        return m_Nodes.back();
    }

    Node2D & NodePool::getNode(const size_t Index)
    {
        assert(Index <= m_Nodes.size());
        return m_Nodes[Index - 1];
    }

    size_t maxNodeIndex()
    {
        NodePool & pool = NodePool::Instance();
        Node2D aNode = *max_element(
          pool.m_Nodes.begin(), pool.m_Nodes.end(), [](const Node2D & a, const Node2D & b) {
              return a.getNodeNumber() < b.getNodeNumber();
          });
        return aNode.getNodeNumber();
    }

    std::vector<double> properties(const Variable t_Property)
    {
        NodePool & pool = NodePool::Instance();   // get the singleton instance
        std::vector<double> aVector;
        for(Node2D & aNode : pool.m_Nodes)
        {
            aVector.push_back(aNode.property(t_Property));
        }
        return aVector;
    }

    void updateNodeValues(const std::vector<double> & t_values,
                          const BaseVariable t_property,
                          bool updatePreviousTimestep)
    {
        NodePool & pool = NodePool::Instance();   // get the singleton instance

#ifdef STL_MULTITHREADING
        std::for_each(std::execution::par_unseq,
                      std::begin(pool.m_Nodes),
                      std::end(pool.m_Nodes),
                      [&](auto && aNode) {
                          const auto nodeNumber = aNode.getNodeNumber() - 1;
                          aNode.setStateProperty(
                            t_property, t_values[nodeNumber], updatePreviousTimestep);
                      });
#else
        for(auto & node : pool.m_Nodes)
        {
            const auto nodeNumber = node.getNodeNumber() - 1;
            node.setStateProperty(t_property, t_values[nodeNumber], updatePreviousTimestep);
        }
#endif
    }

    void NodePool::clear()
    {
        m_Nodes.clear();
    }
}   // namespace HygroThermFEM
