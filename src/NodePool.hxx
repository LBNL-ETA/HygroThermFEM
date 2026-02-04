#pragma once

#include <vector>
#include "Node2D.hxx"

namespace HygroThermFEM
{
    //! \brief Class holds pool of nodes used in one domain.
    //!
    //! Pool of all nodes will be used by different domains that present thermal, mass or pressure
    //! models. NodePool is owned by MultiDomain and must be cleared before new model is used.
    class NodePool
    {
    public:
        NodePool() = default;
        ~NodePool() = default;

        //! \brief Creates node with node index and x, y coordinate and sets starting state
        //! conditions
        //!
        //! \param NodeIndex Node index in finite element domain.
        //! \param x x-coordinate.
        //! \param y y-coordinate
        //! \param state State of variables in the node.
        //! \return Reference to newly created node
        Node2D & createNode(std::size_t NodeIndex,
                            double x,
                            double y,
                            const State & state = State(0, 0, 0, 0));

        //! \brief Returns node at given index.
        //!
        //! \param Index Node index.
        //! \return Reference to node at requested index
        Node2D & getNode(std::size_t Index);

        //! Returns maximum node index from pool of nodes.
        std::size_t maxIndex() const;

        //! Returns state values (temperature, water content or pressure) at all nodes.
        std::vector<double>
          properties(Variable t_Property   //!< Variable for which values are obtained.
          ) const;

        //! Update all node values for given state property.
        //!
        //! @param t_values Vector of new values for given BaseVariable property.
        //! @param t_property Property for which new values will be applied to
        //! @param updatePreviousTimestep Indicates if old value should be stored as previous
        //! timestep value.
        void updateNodeValues(const std::vector<double> & t_values,
                              BaseVariable t_property,
                              bool updatePreviousTimestep = false);

        //! Delete all nodes from NodePool.
        void clear();

    private:
        //! Storage for nodes in NodePool.
        std::vector<Node2D> m_Nodes;
    };

}   // namespace HygroThermFEM
