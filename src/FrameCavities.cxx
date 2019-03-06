#include "FrameCavities.hxx"
#include "MaterialPool.hxx"

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
        auto frameCavities =
          MaterialPool::Instance().getMaterials(MaterialType::FrameCavity_ISO15099);

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

    FrameCavityBoundaries::line::line(const size_t n1, const size_t n2): n1(n1),
                                                               n2(n2) {
    }

    size_t FrameCavityBoundaries::line::getN1() const {
        return n1;
    }

    size_t FrameCavityBoundaries::line::getN2() const {
        return n2;
    }

    bool FrameCavityBoundaries::line::operator<(const line & rhs) const {
        if(n1 < rhs.n1)
            return true;
        if(rhs.n1 < n1)
            return false;
        return n2 < rhs.n2;
    }

    bool FrameCavityBoundaries::line::operator>(const line & rhs) const {
        return rhs < *this;
    }

    bool FrameCavityBoundaries::line::operator<=(const line & rhs) const {
        return !(rhs < *this);
    }

    bool FrameCavityBoundaries::line::operator>=(const line & rhs) const {
        return !(*this < rhs);
    }
}   // namespace HygroThermFEM
