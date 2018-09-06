#include <cassert>
#include <algorithm>

#include "NodePool.hxx"

namespace MoisThermFEM
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
        auto aNode = Node2D(t_NodeNumber, t_x, t_y, t_Prop);
        m_Nodes.push_back(aNode);
        return m_Nodes.back();
    }

    Node2D & NodePool::getNode(const size_t Index)
    {
        assert(Index <= m_Nodes.size());
        return m_Nodes[Index - 1];
    }

    size_t NodePool::maxIndex() const
    {
        Node2D aNode =
          *max_element(m_Nodes.begin(), m_Nodes.end(), [](const Node2D & a, const Node2D & b) {
              return a.getNodeNumber() < b.getNodeNumber();
          });
        return aNode.getNodeNumber();
    }

    std::vector<double> NodePool::nodeProperties(Property t_Property) const
    {
        std::vector<double> aVector;
        for(const Node2D & aNode : m_Nodes)
        {
            aVector.push_back(aNode.getProperty(t_Property));
        }
        return aVector;
    }

    std::vector<double> NodePool::waterContent(WaterContent content) const
    {
        std::vector<double> aVector;
        for(const Node2D & aNode : m_Nodes)
        {
            aVector.push_back(aNode.waterContent(content));
        }
        return aVector;
    }

    void NodePool::updateNodeValues(const std::vector<double> & t_values, const Property t_property)
    {
        assert(m_Nodes.size() == t_values.size());
        for(std::size_t i = 0; i < t_values.size(); ++i)
        {
            m_Nodes[i].setProperty(t_property, t_values[i]);
        }
    }

    void NodePool::clear()
    {
        m_Nodes.clear();
    }

}   // namespace MoisThermFEM