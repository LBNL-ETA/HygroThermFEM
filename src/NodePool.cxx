#include <assert.h>
#include <algorithm>

#include "NodePool.hxx"
#include "Node2D.hxx"

namespace Conrad {

    NodePool &NodePool::Instance() {
        static NodePool m_Instance;
        return m_Instance;
    }

    Node2D &NodePool::createNode(
            size_t const t_NodeNumber,
            double const t_x,
            double const t_y,
            double const t_temperature) {
        auto aNode = Node2D(t_NodeNumber, t_x, t_y, t_temperature);
        m_Nodes.push_back(aNode);
        return m_Nodes.back();
    }

    Node2D &NodePool::getNode(size_t const Index) {
        assert(Index < m_Nodes.size());
        return m_Nodes[Index];
    }

    size_t NodePool::maxIndex() const {
        Node2D aNode = *max_element(m_Nodes.begin(), m_Nodes.end(),
                                    [](const Node2D &a, const Node2D &b) { return a.nodeNumber < b.nodeNumber; });
        return aNode.nodeNumber;
    }

    std::vector<double> NodePool::nodeTemperatures() const {
        std::vector<double> aVector;
        for (const Node2D &aNode : m_Nodes) {
            aVector.push_back(aNode.temperature);
        }
        return aVector;
    }

    void NodePool::clear() {
        m_Nodes.clear();
    }

    NodePool::NodePool() {

    }

    NodePool::~NodePool() {

    }

}