#include <cassert>
#include <algorithm>
#include <ranges>

#include "lbnl/algorithm.hxx"

#include "Nodes.hxx"

#ifdef STL_MULTITHREADING
// par rather than par_unseq: every node update recomputes its water content, which allocates,
// and allocation is not one of the vectorisation-safe operations par_unseq permits. Each node
// writes only to itself, so the result does not depend on how the work is scheduled.
#include <execution>
#endif

namespace HygroThermFEM
{
    Node2D & Nodes::createNode(const NodeParams & params)
    {
        const auto index = params.index == 0 ? m_Nodes.size() + 1 : params.index;
        m_Nodes.emplace_back(index, params.x, params.y, params.state);
        return m_Nodes.back();
    }

    Node2D & Nodes::getNode(const size_t Index)
    {
        assert(Index <= m_Nodes.size());
        return m_Nodes[Index - 1];
    }

    size_t Nodes::maxIndex() const
    {
        // Iterator, not a copy: Node2D owns several node-local containers, so materialising
        // one here allocated on every call -- and this runs twice per matrix assembly.
        const auto maxNode = std::ranges::max_element(
          m_Nodes, {}, [](const Node2D & node) { return node.getNodeNumber(); });

        return maxNode->getNodeNumber();
    }

    std::vector<double> Nodes::properties(const Variable t_Property) const
    {
        return lbnl::transform_to_vector(
          m_Nodes, [&](const Node2D & aNode) { return aNode.property(t_Property); });
    }

    std::vector<double> Nodes::properties(const Variable t_Property,
                                          const Timestep t_Timestep) const
    {
        return lbnl::transform_to_vector(m_Nodes, [&](const Node2D & aNode) {
            return aNode.property(t_Property, t_Timestep);
        });
    }

    void Nodes::updateNodeTemperatures(const std::vector<double> & values,
                                        bool updatePreviousTimestep)
    {
        assert(m_Nodes.size() == values.size());

#ifdef STL_MULTITHREADING
        std::for_each(
          std::execution::par, std::begin(m_Nodes), std::end(m_Nodes), [&](auto && aNode) {
              const auto nodeNumber = aNode.getNodeNumber() - 1;
              aNode.setTemperature(values[nodeNumber], updatePreviousTimestep);
          });
#else
        for(auto & node : m_Nodes)
        {
            const auto nodeNumber = node.getNodeNumber() - 1;
            node.setTemperature(values[nodeNumber], updatePreviousTimestep);
        }
#endif
    }

    void Nodes::updateNodeHumidities(const std::vector<double> & values,
                                      bool updatePreviousTimestep)
    {
        assert(m_Nodes.size() == values.size());

#ifdef STL_MULTITHREADING
        std::for_each(
          std::execution::par, std::begin(m_Nodes), std::end(m_Nodes), [&](auto && aNode) {
              const auto nodeNumber = aNode.getNodeNumber() - 1;
              aNode.setHumidity(values[nodeNumber], updatePreviousTimestep);
          });
#else
        for(auto & node : m_Nodes)
        {
            const auto nodeNumber = node.getNodeNumber() - 1;
            node.setHumidity(values[nodeNumber], updatePreviousTimestep);
        }
#endif
    }

    void Nodes::updateNodeLiquidPercents(const std::vector<double> & values,
                                         bool updatePreviousTimestep)
    {
        assert(m_Nodes.size() == values.size());

#ifdef STL_MULTITHREADING
        std::for_each(
          std::execution::par, std::begin(m_Nodes), std::end(m_Nodes), [&](auto && aNode) {
              const auto nodeNumber = aNode.getNodeNumber() - 1;
              aNode.setLiquidPercent(values[nodeNumber], updatePreviousTimestep);
          });
#else
        for(auto & node : m_Nodes)
        {
            const auto nodeNumber = node.getNodeNumber() - 1;
            node.setLiquidPercent(values[nodeNumber], updatePreviousTimestep);
        }
#endif
    }

    void Nodes::clear()
    {
        m_Nodes.clear();
    }

}   // namespace HygroThermFEM
