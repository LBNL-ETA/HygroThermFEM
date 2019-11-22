#include <KeffCavity.hxx>

#include "FrameCavities.hxx"
#include "MaterialPool.hxx"
#include "NodePool.hxx"
#include "Common.hxx"

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////
    ///  EquivalentFrameCavity
    ///////////////////////////////////////////////////////////////////////////////

    EquivalentFrameCavity::EquivalentFrameCavity(const std::vector<size_t> & nodes,
                                                 IGas & gas,
                                                 FenestrationCommon::GravityVector gravityVector) :
        m_Segments(buildSegments(nodes)),
        m_SideSegments(groupSegmentSides(m_Segments)),
        m_Side{
          {Side::Top, {0, 0}}, {Side::Bottom, {0, 0}}, {Side::Left, {0, 0}}, {Side::Right, {0, 0}}},
        m_Area(calcArea()),
        m_Size(calcSize(m_Area)),
        m_Gas(gas),
        m_GravityVector(std::move(gravityVector))
    {
        calcSideEmissivities();
        updateSideTemperatures();
    }

    void EquivalentFrameCavity::update()
    {
        updateSideTemperatures();
        const auto hfDirection = heatFlowDirection();
        const auto side1 = getSide1(hfDirection);
        const auto side2 = getSide2(hfDirection);
        const auto standard = m_Gas.standard();
        const auto radCalc{KeffCavity::RadiationCalculation::Yes};
        auto thermalConductivity{0.0};
        switch(standard)
        {
            case CavityStandard::ISO15099:
            {
                const auto jambHeight = 1;
                const auto pressure = 101325;
                KeffCavity::CavityISO10599 cavity(hfDirection,
                                                  m_Size.L,
                                                  m_Size.H,
                                                  jambHeight,
                                                  side1,
                                                  side2,
                                                  pressure,
                                                  m_GravityVector,
                                                  radCalc,
                                                  Gases::CGas());
                thermalConductivity = cavity.effectiveConductivity();
                break;
            }
            case CavityStandard::CEN:
            {
                KeffCavity::CavityCEN cavity(
                  hfDirection, m_Size.L, m_Size.H, m_Area, side1, side2, radCalc);
                thermalConductivity = cavity.effectiveConductivity();
                break;
            }
            default:;
        }
        m_Gas.updateThermalConductivity(thermalConductivity);
    }

    std::vector<EquivalentFrameCavity::Segment>
      EquivalentFrameCavity::buildSegments(const std::vector<size_t> & nodes)
    {
        std::vector<Segment> segments;
        for(size_t i = 0u; i < nodes.size(); ++i)
        {
            const auto firstIndex = i == 0 ? nodes.size() - 1 : i - 1;
            const auto & node1 = NodePool::Instance().getNode(nodes[firstIndex]);
            const auto & node2 = NodePool::Instance().getNode(nodes[i]);
            auto emissivity{0.0};
            const auto materialName = findCommonMaterial(node1, node2);
            if(!materialName.empty())
            {
                const auto & material = MaterialPool::Instance().material(materialName);
                emissivity = material.emissivity();
            }
            segments.emplace_back(node1, node2, emissivity);
        }
        return segments;
    }

    KeffCavity::CavitySide EquivalentFrameCavity::getSide1(KeffCavity::ScreenFlow screenFlow)
    {
        std::map<KeffCavity::ScreenFlow, Side> flowMap{
          {KeffCavity::ScreenFlow::Left, Side::Right},
          {KeffCavity::ScreenFlow::Right, Side::Left},
          {KeffCavity::ScreenFlow::Up, Side::Bottom},
          {KeffCavity::ScreenFlow::Down, Side::Top},
        };

        return m_Side.at(flowMap.at(screenFlow));
    }

    KeffCavity::CavitySide EquivalentFrameCavity::getSide2(KeffCavity::ScreenFlow screenFlow)
    {
        std::map<KeffCavity::ScreenFlow, Side> flowMap{
          {KeffCavity::ScreenFlow::Left, Side::Left},
          {KeffCavity::ScreenFlow::Right, Side::Right},
          {KeffCavity::ScreenFlow::Up, Side::Top},
          {KeffCavity::ScreenFlow::Down, Side::Bottom},
        };

        return m_Side.at(flowMap.at(screenFlow));
    }

    std::string EquivalentFrameCavity::findCommonMaterial(const Node2D & node1,
                                                          const Node2D & node2)
    {
        std::string name;
        auto node1Materials = node1.getSolidMaterialNames();
        auto node2Materials = node2.getSolidMaterialNames();

        for(const auto & mat1 : node1Materials)
        {
            for(const auto & mat2 : node2Materials)
            {
                if(mat1 == mat2)
                {
                    name = mat1;
                    break;
                }
            }
        }

        return name;
    }

    double EquivalentFrameCavity::calcArea() const
    {
        double area{0};
        for(const auto & segment : m_Segments)
        {
            area += segment.crossCalc();
        }
        return 0.5 * area;
    }

    EquivalentFrameCavity::Size EquivalentFrameCavity::calcSize(const double area) const
    {
        double maxX = m_Segments[0].firstNode().X();
        double minX = m_Segments[0].firstNode().X();
        double maxY = m_Segments[0].firstNode().Y();
        double minY = m_Segments[0].firstNode().Y();

        for(size_t i = 1u; i < m_Segments.size(); ++i)
        {
            maxX = std::max(maxX, m_Segments[i].firstNode().X());
            minX = std::min(minX, m_Segments[i].firstNode().X());
            maxY = std::max(maxY, m_Segments[i].firstNode().Y());
            minY = std::min(minY, m_Segments[i].firstNode().Y());
        }

        const auto XSize = maxX - minX;
        const auto YSize = maxY - minY;
        const auto ratio = XSize / YSize;

        const auto H = std::sqrt(area / ratio);
        const auto L = ratio * H;

        return {L, H};
    }

    std::map<EquivalentFrameCavity::Side, std::vector<EquivalentFrameCavity::Segment>>
      EquivalentFrameCavity::groupSegmentSides(const std::vector<Segment> & segments)
    {
        std::map<Side, std::vector<Segment>> result{{Side::Top, std::vector<Segment>()},
                                                    {Side::Bottom, std::vector<Segment>()},
                                                    {Side::Left, std::vector<Segment>()},
                                                    {Side::Right, std::vector<Segment>()}};
        for(const auto & segment : segments)
        {
            result.at(segment.side()).push_back(segment);
        }
        return result;
    }

    void EquivalentFrameCavity::calcSideEmissivities()
    {
        std::map<Side, double> length{
          {Side::Top, 0}, {Side::Bottom, 0}, {Side::Left, 0}, {Side::Right, 0}};
        std::map<Side, double> emissLength{
          {Side::Top, 0}, {Side::Bottom, 0}, {Side::Left, 0}, {Side::Right, 0}};
        for(const auto & side : m_SideSegments)
        {
            const auto aSide = side.first;
            for(const auto & segment : side.second)
            {
                length.at(aSide) += segment.length();
                emissLength.at(aSide) += segment.length() * segment.emissivity();
            }
        }

        for(const auto side : EnumSide())
        {
            if(length.at(side) != 0)
            {
                m_Side.at(side).emissivity = emissLength.at(side) / length.at(side);
            }
        }
    }

    void EquivalentFrameCavity::updateSideTemperatures()
    {
        std::map<Side, double> length{
          {Side::Top, 0}, {Side::Bottom, 0}, {Side::Left, 0}, {Side::Right, 0}};
        std::map<Side, double> temperatureLength{
          {Side::Top, 0}, {Side::Bottom, 0}, {Side::Left, 0}, {Side::Right, 0}};
        for(const auto & side : m_SideSegments)
        {
            const auto aSide = side.first;
            for(const auto & segment : side.second)
            {
                length.at(aSide) += segment.length();
                temperatureLength.at(aSide) += segment.length() * toKelvin(segment.temperature());
            }
        }

        // TODO: Standard enumerator won't work. Fix this later.
        for(const auto side : EnumSide())
        {
            if(length.at(side) != 0)
            {
                m_Side.at(side).temperature = temperatureLength.at(side) / length.at(side);
            }
        }
    }

    KeffCavity::ScreenFlow EquivalentFrameCavity::heatFlowDirection() const
    {
        const auto leftTemp = m_Side.at(Side::Left).temperature;
        const auto rightTemp = m_Side.at(Side::Right).temperature;
        const auto topTemp = m_Side.at(Side::Top).temperature;
        const auto bottomTemp = m_Side.at(Side::Bottom).temperature;
        KeffCavity::ScreenFlow result{KeffCavity::ScreenFlow::Right};
        if(std::abs(topTemp - bottomTemp) > std::abs(leftTemp - rightTemp))
        {
            if(topTemp > bottomTemp)
            {
                result = KeffCavity::ScreenFlow::Down;
            }
            else
            {
                result = KeffCavity::ScreenFlow::Up;
            }
        }
        if(rightTemp > leftTemp)
        {
            result = KeffCavity::ScreenFlow::Left;
        }

        return result;
    }

    EquivalentFrameCavity::Segment::Segment(const Node2D & node1,
                                            const Node2D & node2,
                                            double emissivity) :
        node1(node1),
        node2(node2),
        m_Emissivity(emissivity),
        m_Length(
          std::sqrt(std::pow(node1.X() - node2.X(), 2) + std::pow(node1.Y() - node2.Y(), 2))),
        m_Side(calcSide(node1, node2))
    {}

    EquivalentFrameCavity::Side EquivalentFrameCavity::Segment::calcSide(const Node2D & n1,
                                                                         const Node2D & n2) const
    {
        const auto PI = 4 * std::atan(1);
        const auto angle = std::atan2((n2.Y() - n1.Y()), (n2.X() - n1.X()));
        Side aSide{Side::Left};
        const auto sectorAngle{PI / 4};
        // TODO: This is processing for counter clock wise. Need mechanism that will catch clockwise
        // direction as well
        if(angle < sectorAngle && angle > -sectorAngle)
        {
            aSide = Side::Bottom;
        }
        else if(angle >= sectorAngle && angle < 3 * sectorAngle)
        {
            aSide = Side::Right;
        }
        else if(angle < -sectorAngle && angle > -3 * sectorAngle)
        {
            aSide = Side::Left;
        }
        else if((angle > 3 * sectorAngle || angle <= -3 * sectorAngle))
        {
            aSide = Side::Top;
        }
        return aSide;
    }

    double EquivalentFrameCavity::Segment::emissivity() const
    {
        return m_Emissivity;
    }

    double EquivalentFrameCavity::Segment::length() const
    {
        return m_Length;
    }

    double EquivalentFrameCavity::Segment::temperature() const
    {
        return 0.5
               * (node1.property(Variable::temperature) + node2.property(Variable::temperature));
    }

    double EquivalentFrameCavity::Segment::crossCalc() const
    {
        return node1.X() * node2.Y() - node1.Y() * node2.X();
    }

    const Node2D & EquivalentFrameCavity::Segment::firstNode() const
    {
        return node1;
    }

    const EquivalentFrameCavity::Side & EquivalentFrameCavity::Segment::side() const
    {
        return m_Side;
    }

    ///////////////////////////////////////////////////////////////////////////////
    ///  EquivalentFrameCavities
    ///////////////////////////////////////////////////////////////////////////////

    EquivalentFrameCavities::EquivalentFrameCavities(const ElementsLinear2D & m_Elements) :
        m_Elements(m_Elements)
    {
        createEquivalentFrameCavities();
    }

    void EquivalentFrameCavities::createEquivalentFrameCavities()
    {
        auto frameCavities = MaterialPool::Instance().getGases();

        for(const auto & frameCavity : frameCavities)
        {
            std::vector<std::vector<size_t>> elementNodes;
            for(const auto & element : m_Elements.elements())
            {
                if(element->getMaterial().name() == frameCavity)
                {
                    elementNodes.push_back(element->nodeIndexes());
                }
            }
            auto edges = getEdges(elementNodes);
            auto & gas = MaterialPool::Instance().gas(frameCavity);
            m_Cavities.emplace_back(edgeNodesOrdered(edges), gas);
        }
    }

    std::set<EquivalentFrameCavities::line>
      EquivalentFrameCavities::getEdges(const std::vector<std::vector<size_t>> & elNodes)
    {
        std::map<std::set<size_t>, size_t> edges;
        std::set<line> allEdges;

        for(auto & n : elNodes)
        {
            for(size_t i = 0u; i < n.size(); ++i)
            {
                const auto index1 = n[i];
                auto index2 = n[0u];
                if(i != n.size() - 1)
                {
                    index2 = n[i + 1u];
                }
                std::set<size_t> edge({index1, index2});
                if(edges.find(edge) != edges.end())
                {
                    edges[edge]++;
                    allEdges.erase(line(index1, index2));
                    allEdges.erase(line(index2, index1));
                }
                else
                {
                    edges[edge] = 1u;
                    allEdges.insert(line(index1, index2));
                }
            }
        }

        return allEdges;
    }

    std::vector<size_t> EquivalentFrameCavities::edgeNodesOrdered(std::set<line> & allEdges)
    {
        std::vector<size_t> boundaryLine;

        line first = *allEdges.begin();
        boundaryLine.push_back(first.getN1());
        allEdges.erase(first);

        while(!allEdges.empty())
        {
            auto second = std::find_if(allEdges.begin(), allEdges.end(), [&first](const line & l) {
                return first.getN2() == l.getN1();
            });
            first = *second;
            boundaryLine.push_back(first.getN1());
            allEdges.erase(first);
        }

        return boundaryLine;
    }

    void EquivalentFrameCavities::update()
    {
        for(auto & cavity : m_Cavities)
        {
            cavity.update();
        }
    }

    EquivalentFrameCavities::line::line(const size_t n1, const size_t n2) : n1(n1), n2(n2)
    {}

    size_t EquivalentFrameCavities::line::getN1() const
    {
        return n1;
    }

    size_t EquivalentFrameCavities::line::getN2() const
    {
        return n2;
    }

    bool EquivalentFrameCavities::line::operator<(const line & rhs) const
    {
        if(n1 < rhs.n1)
            return true;
        if(rhs.n1 < n1)
            return false;
        return n2 < rhs.n2;
    }

    bool EquivalentFrameCavities::line::operator>(const line & rhs) const
    {
        return rhs < *this;
    }

    bool EquivalentFrameCavities::line::operator<=(const line & rhs) const
    {
        return !(rhs < *this);
    }

    bool EquivalentFrameCavities::line::operator>=(const line & rhs) const
    {
        return !(*this < rhs);
    }
}   // namespace HygroThermFEM
