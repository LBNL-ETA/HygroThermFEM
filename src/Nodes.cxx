#include <cassert>
#include <algorithm>

#include "Nodes.hxx"

#ifdef STL_MULTITHREADING
#include <execution>
#endif

namespace HygroThermFEM
{
    Node2D & Nodes::createNode(const std::size_t t_NodeNumber,
                                  const double t_x,
                                  const double t_y,
                                  const State & t_Prop)
    {
        m_Nodes.emplace_back(t_NodeNumber, t_x, t_y, t_Prop);
        return m_Nodes.back();
    }

    Node2D & Nodes::createNode(const double x, const double y, const State & state)
    {
        const auto nodeIndex = m_Nodes.size() + 1;
        return createNode(nodeIndex, x, y, state);
    }

    Node2D & Nodes::getNode(const size_t Index)
    {
        assert(Index <= m_Nodes.size());
        return m_Nodes[Index - 1];
    }

    size_t Nodes::maxIndex() const
    {
        Node2D aNode =
          *max_element(m_Nodes.begin(), m_Nodes.end(), [](const Node2D & a, const Node2D & b) {
              return a.getNodeNumber() < b.getNodeNumber();
          });
        return aNode.getNodeNumber();
    }

    std::vector<double> Nodes::properties(const Variable t_Property) const
    {
        std::vector<double> aVector;
        for(const Node2D & aNode : m_Nodes)
        {
            aVector.push_back(aNode.property(t_Property));
        }
        return aVector;
    }

    void Nodes::updateNodeValues(const std::vector<double> & t_values,
                                    const BaseVariable t_property,
                                    bool updatePreviousTimestep)
    {
        assert(m_Nodes.size() == t_values.size());

#ifdef STL_MULTITHREADING
        std::for_each(
          std::execution::par_unseq, std::begin(m_Nodes), std::end(m_Nodes), [&](auto && aNode) {
              const auto nodeNumber = aNode.getNodeNumber() - 1;
              aNode.setStateProperty(t_property, t_values[nodeNumber], updatePreviousTimestep);
          });
#else
        for(auto & node: m_Nodes)
        {
            const auto nodeNumber = node.getNodeNumber() - 1;
            node.setStateProperty(t_property, t_values[nodeNumber], updatePreviousTimestep);
        }
#endif
    }

    void Nodes::clear()
    {
        m_Nodes.clear();
    }

}   // namespace HygroThermFEM
