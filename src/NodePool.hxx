#pragma once

#include <vector>
#include "Node2D.hxx"

namespace MoisThermFEM
{
    class NodePool
    {
    public:
        static NodePool & Instance();

        Node2D & createNode(const std::size_t t_NodeNumber,
                            const double t_x,
                            const double t_y,
                            const State & t_Prop = State(0, 0, 0, 0));

        Node2D & getNode(std::size_t const Index);

        std::size_t maxIndex() const;
        std::vector<double> nodeProperties(Property t_Property) const;
        std::vector<double> waterContent() const;
        std::vector<double> liquidContent() const;
        std::vector<double> vaporContent() const;
        std::vector<double> iceContent() const;

        void updateNodeValues(const std::vector<double> & t_values, const Property t_property);

        void clear();

    private:
        NodePool() = default;
        ~NodePool() = default;

        std::vector<Node2D> m_Nodes;
    };

}   // namespace MoisThermFEM