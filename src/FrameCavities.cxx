#include "FrameCavities.hxx"
#include "MaterialPool.hxx"
#include "NodePool.hxx"

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////
    ///  FrameCavity
    ///////////////////////////////////////////////////////////////////////////////

    FrameCavity::FrameCavity(const double effectiveConductivity,
                             const double horizontalDimension,
                             const double verticalDimension,
                             const double nusseltNumber,
                             const HygroThermFEM::FrameCavityScreenHeatFlow heatFlowDirection,
                             const double emissivity1,
                             const double temperature1,
                             const double emissivity2,
                             const double temperature2) :
        m_EffectiveConductivity(effectiveConductivity),
        m_HorizontalDimension(horizontalDimension),
        m_VerticalDimension(verticalDimension),
        m_NusseltNumber(nusseltNumber),
        m_HeatFlowDirection(heatFlowDirection),
        m_Side1{temperature1, emissivity1},
        m_Side2{temperature2, emissivity2}
    {}

    ///////////////////////////////////////////////////////////////////////////////
    ///  FrameCavityBoundaries
    ///////////////////////////////////////////////////////////////////////////////

    FrameCavityBoundaries::FrameCavityBoundaries(const ElementsLinear2D & m_Elements) :
        m_Elements(m_Elements)
    {
        calculateEquivalentFrameCavities();
    }

    void FrameCavityBoundaries::calculateEquivalentFrameCavities()
    {
        auto frameCavities = MaterialPool::Instance().getMaterials(MaterialType::Gas);

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
            m_BoundaryNodes[frameCavity] = edgeNodesOrdered(edges);
        }
    }

    std::set<FrameCavityBoundaries::line>
      FrameCavityBoundaries::getEdges(const std::vector<std::vector<size_t>> & elNodes)
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

    std::vector<size_t> FrameCavityBoundaries::edgeNodesOrdered(std::set<line> & allEdges)
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

    const std::vector<size_t> &
      FrameCavityBoundaries::boundaryNodes(const std::string & frameCavityName) const
    {
        return m_BoundaryNodes.at(frameCavityName);
    }

    FrameCavityBoundaries::line::line(const size_t n1, const size_t n2) : n1(n1), n2(n2)
    {}

    size_t FrameCavityBoundaries::line::getN1() const
    {
        return n1;
    }

    size_t FrameCavityBoundaries::line::getN2() const
    {
        return n2;
    }

    bool FrameCavityBoundaries::line::operator<(const line & rhs) const
    {
        if(n1 < rhs.n1)
            return true;
        if(rhs.n1 < n1)
            return false;
        return n2 < rhs.n2;
    }

    bool FrameCavityBoundaries::line::operator>(const line & rhs) const
    {
        return rhs < *this;
    }

    bool FrameCavityBoundaries::line::operator<=(const line & rhs) const
    {
        return !(rhs < *this);
    }

    bool FrameCavityBoundaries::line::operator>=(const line & rhs) const
    {
        return !(*this < rhs);
    }

    RectangularizedCavity::RectangularizedCavity(const std::vector<size_t> & nodes) :
        m_Segments(buildSegments(nodes)),
        m_Area(area())
    {}

    std::vector<RectangularizedCavity::Segment>
      RectangularizedCavity::buildSegments(const std::vector<size_t> & nodes)
    {
        std::vector<Segment> segments;
        for(size_t i = 1u; i < nodes.size(); ++i)
        {
            const auto & node1 = NodePool::Instance().getNode(i - 1);
            const auto & node2 = NodePool::Instance().getNode(i);
            const auto materialName = findCommonMaterial(node1, node2);
            const auto & material = MaterialPool::Instance().material(materialName);
            segments.emplace_back(node1, node2, material.emissivity());
        }
        return segments;
    }

    std::string RectangularizedCavity::findCommonMaterial(const Node2D & node1,
                                                          const Node2D & node2) const
    {
        std::string name;
        auto node1Materials = node1.getMaterialNames(MaterialType::Solid);
        auto node2Materials = node2.getMaterialNames(MaterialType::Solid);

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

    double RectangularizedCavity::area() const
    {
        double area{0};
        for(const auto & segment : m_Segments)
        {
            area += segment.crossCalc();
        }
        return 0.5 * area;
    }

    RectangularizedCavity::Segment::Segment(const Node2D & node1,
                                            const Node2D & node2,
                                            double emissivity) :
        node1(node1),
        node2(node2),
        m_Emissivity(emissivity),
        m_Length(
          std::sqrt(std::pow(node1.X() - node2.X(), 2) + std::pow(node1.Y() - node2.Y(), 2))),
        m_Side(calcSide(node1, node2))
    {}

    RectangularizedCavity::Side RectangularizedCavity::Segment::calcSide(const Node2D & node1,
                                                                         const Node2D & node2) const
    {
        auto angle =
          node1.X() != node2.X() ? atan(node2.Y() - node1.Y() / (node2.X() - node1.X())) : M_PI / 2;
        HygroThermFEM::RectangularizedCavity::Side aSide{Side::Left};
        const auto sectorAngle{M_PI / 4};
        if(angle >= sectorAngle && angle < 3 * sectorAngle)
        {
            aSide = HygroThermFEM::RectangularizedCavity::Side::Left;
        }
        else if(angle >= 3 * sectorAngle && angle < 5 * sectorAngle)
        {
            aSide = HygroThermFEM::RectangularizedCavity::Side::Bottom;
        }
        else if(angle >= 5 * sectorAngle && angle < 7 * sectorAngle)
        {
            aSide = HygroThermFEM::RectangularizedCavity::Side::Right;
        }
        else if((angle < sectorAngle && angle >= 7 * sectorAngle))
        {
            aSide = HygroThermFEM::RectangularizedCavity::Side::Top;
        }
        return aSide;
    }

    double RectangularizedCavity::Segment::emissivity() const
    {
        return m_Emissivity;
    }

    double RectangularizedCavity::Segment::node1Temperature() const
    {
        return node1.property(Variable::temperature);
    }

    double RectangularizedCavity::Segment::node2Temperature() const
    {
        return node2.property(Variable::temperature);
    }

    double RectangularizedCavity::Segment::length() const
    {
        return m_Length;
    }

    double RectangularizedCavity::Segment::averageTemperature() const
    {
        return 0.5
               * (node1.property(Variable::temperature) + node2.property(Variable::temperature));
    }

    double RectangularizedCavity::Segment::crossCalc() const
    {
        return node1.X() * node2.Y() - node1.Y() * node2.X();
    }
}   // namespace HygroThermFEM
