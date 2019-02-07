#pragma once

#include <vector>
#include "Node2D.hxx"

namespace HygroThermFEM
{
    //! \brief Singleton class holds poll of nodes used in one domain.
    //!
    //! Pool of all nodes will be used by different domains that present thermal, mass or pressure
    //! models. NodePool must be cleared up before new model is used.
    class NodePool
    {
    public:
        //! Access to NodePool singleton
        static NodePool & Instance();

        //! Creates node with node index and x, y coordinate and sets starting state conditions
        Node2D &
          createNode(std::size_t NodeIndex,   //!< Node index in finite element domain.
                     double x,                //!< x-coordinate.
                     double y,                //!< y-coordinate
                     const State & state = State(0, 0, 0, 0)   //!< State of variables in the node.
          );

        //! Returns node at given index.
        Node2D & getNode(std::size_t Index   //!< Node index.
        );

        //! Returns maximum node index from pool of nodes.
        std::size_t maxIndex() const;

        //! Returns state values (temperature, water content or pressure) at all nodes.
        std::vector<double>
          properties(Variable t_Property   //!< Variable for which values are obtained.
          );

        //! Update all node values for given state property.
        void updateNodeValues(
          const std::vector<double> & t_values,   //!< Vector of new values for given BaseVariable.
          BaseVariable t_property   //!< BaseVariable for which new values will be applied to.
        );

        //! Delete all nodes from NodePool.
        void clear();

    private:
        NodePool() = default;
        ~NodePool() = default;

        //! Storage for nodes in NodePool.
        std::vector<Node2D> m_Nodes;
    };

}   // namespace HygroThermFEM
