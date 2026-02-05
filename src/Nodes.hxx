#pragma once

#include <vector>
#include "Node2D.hxx"

namespace HygroThermFEM
{
    //! \brief Parameters for creating a node with C++20 designated initializers.
    struct NodeParams
    {
        std::size_t index = 0;
        double x = 0.0;
        double y = 0.0;
        State state = {};
    };

    //! \brief Class holds pool of nodes used in one domain.
    //!
    //! Pool of all nodes will be used by different domains that present thermal, mass or pressure
    //! models. NodePool is owned by MultiDomain and must be cleared before new model is used.
    class Nodes
    {
    public:
        Nodes() = default;
        ~Nodes() = default;

        //! \brief Creates node using C++20 designated initializers.
        //!
        //! Example: createNode({.index = 1, .x = 0.0, .y = 0.5, .state = {.temperature = 20}})
        //! \param params Node parameters.
        //! \return Reference to newly created node
        Node2D & createNode(const NodeParams & params);

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

        //! Update all node temperature values.
        //!
        //! @param values Vector of new temperature values.
        //! @param updatePreviousTimestep Indicates if old value should be stored as previous
        //! timestep value.
        void updateNodeTemperatures(const std::vector<double> & values,
                                    bool updatePreviousTimestep = false);

        //! Update all node humidity values.
        //!
        //! @param values Vector of new humidity values.
        //! @param updatePreviousTimestep Indicates if old value should be stored as previous
        //! timestep value.
        void updateNodeHumidities(const std::vector<double> & values,
                                  bool updatePreviousTimestep = false);

        //! Delete all nodes from NodePool.
        void clear();

    private:
        //! Storage for nodes in NodePool.
        std::vector<Node2D> m_Nodes;
    };

}   // namespace HygroThermFEM
