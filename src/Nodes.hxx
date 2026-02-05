#pragma once

#include <vector>
#include "Node2D.hxx"

namespace HygroThermFEM
{
    //! \brief Class holds pool of nodes used in one domain.
    //!
    //! Pool of all nodes will be used by different domains that present thermal, mass or pressure
    //! models. NodePool is owned by MultiDomain and must be cleared before new model is used.
    class Nodes
    {
    public:
        Nodes() = default;
        ~Nodes() = default;

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
                            const State & state = {.temperature = 0, .humidity = 0, .pressure = 0, .liquidPercent = 0});

        //! \brief Creates node with auto-generated index
        //!
        //! Node index is automatically assigned as the next sequential value.
        //! \param x x-coordinate.
        //! \param y y-coordinate
        //! \param state State of variables in the node.
        //! \return Reference to newly created node
        Node2D & createNode(double x, double y, const State & state = {.temperature = 0, .humidity = 0, .pressure = 0, .liquidPercent = 0});

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
