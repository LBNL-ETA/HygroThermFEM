#include <stdexcept>

#include "Node2D.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint1D
    ////////////////////////////////////////////////////////////////////////////

    LocalPoint1D::LocalPoint1D(double const t_ksi) :
        ksi(t_ksi)
    {}

    LocalPoint1D::LocalPoint1D(LocalPoint1D const & t_LocalPoint)
    {
        ksi = t_LocalPoint.ksi;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///   LocalPoint2D
    ////////////////////////////////////////////////////////////////////////////

    LocalPoint2D::LocalPoint2D(double const t_ksi, double const t_eta) :
        ksi(t_ksi),
        eta(t_eta)
    {}

    LocalPoint2D::LocalPoint2D(LocalPoint2D const & t_LocalPoint)
    {
        ksi = t_LocalPoint.ksi;
        eta = t_LocalPoint.eta;
    }

    INode2D::INode2D(std::size_t t_NodeNumber, double t_x, double t_y) :
        m_NodeNumber(t_NodeNumber),
        m_x(t_x),
        m_y(t_y)
    {}

    bool operator==(const INode2D & first, const INode2D & second)
    {
        bool identical = true;
        identical = identical && first.m_NodeNumber == second.m_NodeNumber;
        identical = identical && first.m_x == second.m_x;
        identical = identical && first.m_y == second.m_y;
        return identical;
    }

    bool operator!=(const INode2D & first, const INode2D & second)
    {
        return !operator==(first, second);
    }

    size_t INode2D::getNodeNumber() const
    {
        return m_NodeNumber;
    }

    double INode2D::X() const
    {
        return m_x;
    }

    double INode2D::Y() const
    {
        return m_y;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///   Node2D
    ////////////////////////////////////////////////////////////////////////////

    Node2D::Node2D(const std::size_t t_NodeNumber,
                   const double t_x,
                   const double t_y,
                   const State & t_State) :
        INode2D(t_NodeNumber, t_x, t_y),
        m_State{{Timestep::Current, t_State}, {Timestep::Previous, t_State}},
        m_Water(calcWaterContent())
    {}

    double Node2D::property(const Variable property, const Timestep iteration) const
    {
        switch(property)
        {
            case Variable::temperature:
                return m_State.at(iteration).getValue(BaseVariable::temperature);
            case Variable::humidity:
                return m_State.at(iteration).getValue(BaseVariable::humidity);
            case Variable::pressure:
                return m_State.at(iteration).getValue(BaseVariable::pressure);
            case Variable::liquidPercent:
                return m_State.at(iteration).getValue(BaseVariable::liquidPercent);
            case Variable::water:
                return waterContent(WaterContent::Water);
            case Variable::liquid:
                return waterContent(WaterContent::Liquid);
            case Variable::vapor:
                return waterContent(WaterContent::Vapor);
            case Variable::ice:
                return waterContent(WaterContent::Ice);
        }
        return 0;
    }

    void Node2D::setStateProperty(const BaseVariable t_Property,
                                  double t_value,
                                  bool updatePreviousValue)
    {
        // First store current to previous iteration and then store current.
        if(updatePreviousValue)
        {
            m_State.at(Timestep::Previous)
                   .setValue(t_Property, m_State.at(Timestep::Current).getValue(t_Property));
        }
        m_State.at(Timestep::Current).setValue(t_Property, t_value);
        updateWaterContent();
    }

    void Node2D::assignMaterial(const IMaterial & material, const double weightingCoefficient)
    {
        const auto & materialName = material.name();
        const auto isAlreadyThere{m_Materials.find(materialName)};
        if(isAlreadyThere != m_Materials.end())
        {
            m_Materials.at(materialName).addWeight(weightingCoefficient);
        }
        m_Materials.emplace(materialName, MaterialContainer{weightingCoefficient, material});
        updateWaterContent();
    }

    double Node2D::waterContent(const WaterContent content) const
    {
        return m_Water.content(content);
    }

    Water Node2D::calcWaterContent()
    {
        Water sum;
        double weighting = 0;
        Water result;
        // Node can end up in several different elements and elements can have different materials.
        // This part of code will check influence (weighting) of different materials on current
        // node. In this way, program will estimate water content of node containing different
        // materials
        for(auto & val : m_Materials)
        {
            auto & material{val.second.material};
            // Since material properties are optional, engine will not be able to update water
            // content in case some properties are missing.
            if(material.hasSorptionCurve() && material.hasPorosity())
            {
                sum += val.second.material.waterContent(*this) * val.second.weightingFactor;
                weighting += val.second.weightingFactor;
            }
        }
        if(weighting != 0)
        {
            result = sum / weighting;
        }
        return result;
    }

    void Node2D::updateWaterContent()
    {
        m_Water = calcWaterContent();
    }

    std::vector<std::string> Node2D::getSolidMaterialNames() const
    {
        std::vector<std::string> names;
        for(const auto & material : m_Materials)
        {
            names.push_back(material.first);
        }
        return names;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///   INodesStorage
    ////////////////////////////////////////////////////////////////////////////

    INodes::INodes(std::initializer_list<std::reference_wrapper<INode2D>> t_Nodes) :
        m_Nodes(t_Nodes)
    {}

    INodes::INodes(Node2D & node1, Node2D & node2)
    {
        m_Nodes.emplace_back(node1);
        m_Nodes.emplace_back(node2);
    }

    INodes::INodes(Node2D & node1, Node2D & node2, Node2D & node3, Node2D & node4)
    {
        m_Nodes.emplace_back(node1);
        m_Nodes.emplace_back(node2);
        m_Nodes.emplace_back(node3);
        m_Nodes.emplace_back(node4);
    }

    std::vector<std::size_t> INodes::getNodeIndexes() const
    {
        std::vector<std::size_t> indexes;
        for(const auto & aNode : m_Nodes)
        {
            indexes.push_back(aNode.get().getNodeNumber());
        }
        return indexes;
    }

    INode2D & INodes::operator[](const std::size_t index) const
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

    std::vector<double> INodes::properties(const Variable property) const
    {
        std::vector<double> result;
        for(const auto & node : m_Nodes)
        {
            result.push_back(node.get().property(property));
        }
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////
    ///   LineNodes2D
    ////////////////////////////////////////////////////////////////////////////

    LineNodes2D::LineNodes2D(Node2D & t_Node1, Node2D & t_Node2) :
        INodes{t_Node1, t_Node2}
    {}

    ////////////////////////////////////////////////////////////////////////////
    ///   QuadrilateralNodes2D
    ////////////////////////////////////////////////////////////////////////////

    QuadrilateralNodes2D::QuadrilateralNodes2D(Node2D & t_Node1,
                                               Node2D & t_Node2,
                                               Node2D & t_Node3,
                                               Node2D & t_Node4) :
        INodes({t_Node1, t_Node2, t_Node3, t_Node4})
    {}
} // namespace HygroThermFEM
