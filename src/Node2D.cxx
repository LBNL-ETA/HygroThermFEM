#include <cassert>

#include "Node2D.hxx"
#include "MaterialPool.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint1D
    ////////////////////////////////////////////////////////////////////////////

    LocalPoint1D::LocalPoint1D(double const t_ksi) : ksi(t_ksi)
    {}

    LocalPoint1D::LocalPoint1D(LocalPoint1D const & t_LocalPoint)
    {
        ksi = t_LocalPoint.ksi;
    }

    ////////////////////////////////////////////////////////////////////////////
    //   LocalPoint2D
    ////////////////////////////////////////////////////////////////////////////

    LocalPoint2D::LocalPoint2D(double const t_ksi, double const t_eta) : ksi(t_ksi), eta(t_eta)
    {}

    LocalPoint2D::LocalPoint2D(LocalPoint2D const & t_LocalPoint)
    {
        ksi = t_LocalPoint.ksi;
        eta = t_LocalPoint.eta;
    }

    ////////////////////////////////////////////////////////////////////////////
    //   Node2D
    ////////////////////////////////////////////////////////////////////////////

    Node2D::Node2D(const std::size_t t_NodeNumber,
                   const double t_x,
                   const double t_y,
                   const State & t_State) :
        m_NodeNumber(t_NodeNumber),
        m_x(t_x),
        m_y(t_y),
        m_State(t_State)
    {}

    bool operator==(const Node2D & first, const Node2D & second)
    {
        bool identical = true;
        identical = identical && first.m_NodeNumber == second.m_NodeNumber;
        identical = identical && first.m_x == second.m_x;
        identical = identical && first.m_y == second.m_y;
        return identical;
    }

    bool operator!=(const Node2D & first, const Node2D & second)
    {
        return !operator==(first, second);
    }

    size_t Node2D::getNodeNumber() const
    {
        return m_NodeNumber;
    }

    double Node2D::X() const
    {
        return m_x;
    }

    double Node2D::Y() const
    {
        return m_y;
    }

    double Node2D::getProperty(const Property t_Property, const Iteration t_Iteration) const
    {
        return m_State.getValue(t_Property, t_Iteration);
    }

    void Node2D::setProperty(const Property t_Property, double t_value)
    {
        m_State.setValue(t_Property, t_value);
    }

    double Node2D::getDeltaProperty(const Property t_Property) const
    {
        return m_State.getDeltaValue(t_Property);
    }

    const State & Node2D::getState() const
    {
        return m_State;
    }

    void Node2D::assignMaterial(const std::string & t_Material, double weightingCoefficient)
    {
        auto & material = MaterialPool::Instance().material(t_Material);
        m_Materials.emplace(weightingCoefficient, material);
    }

    double Node2D::waterContent(WaterContent content) const
    {
        double sum = 0.0;
        double weighting = 0;
        for(auto & val : m_Materials)
        {
            sum += val.second.get().waterContent(m_State, content) * val.first;
            weighting += val.first;
        }
        return sum / weighting;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///   INodesStorage
    ////////////////////////////////////////////////////////////////////////////

    INodes::INodes(std::initializer_list<Node2D> t_Nodes) : m_Nodes(t_Nodes)
    {}

    Node2D & INodes::getNode(const std::size_t Index)
    {
        assert(Index < m_Nodes.size());
        return m_Nodes[Index];
    }

    std::vector<std::size_t> INodes::getNodeIndexes() const
    {
        std::vector<std::size_t> indexes;
        for(const auto & aNode : m_Nodes)
        {
            indexes.push_back(aNode.getNodeNumber());
        }
        return indexes;
    }

    Node2D INodes::operator[](const std::size_t index) const
    {
        if(index >= m_Nodes.size())
        {
            throw std::overflow_error("Index is higher than number of nodes.");
        }
        return m_Nodes[index];
    }

    std::size_t INodes::size() const
    {
        return m_Nodes.size();
    }

	Node2D & INodes::operator[](const size_t index)
	{
		return m_Nodes[index];
	}

	////////////////////////////////////////////////////////////////////////////
    ///   LineNodes2D
    ////////////////////////////////////////////////////////////////////////////

    LineNodes2D::LineNodes2D(const Node2D & t_Node1, const Node2D & t_Node2) :
        INodes{t_Node1, t_Node2}
    {}

    ////////////////////////////////////////////////////////////////////////////
    ///   QuadrilateralNodes2D
    ////////////////////////////////////////////////////////////////////////////

    QuadrilateralNodes2D::QuadrilateralNodes2D(const Node2D & t_Node1,
                                               const Node2D & t_Node2,
                                               const Node2D & t_Node3,
                                               const Node2D & t_Node4) :
        INodes({t_Node1, t_Node2, t_Node3, t_Node4})
    {}
}   // namespace MoisThermFEM
