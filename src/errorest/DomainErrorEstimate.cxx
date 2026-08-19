#include <map>
#include <string>

#include "DomainErrorEstimate.hxx"

#include "Domain.hxx"
#include "Element2D.hxx"
#include "Elements2D.hxx"
#include "Material.hxx"
#include "Node2D.hxx"

namespace HygroThermFEM::errorest
{
    namespace
    {
        //! Stable per-material id so the recovery segregates patches by material, mirroring the
        //! legacy esterror (which recovers per material / vertex-use). Region topology is
        //! all-or-nothing in the estimator, so every element gets a subdomain.
        class MaterialGrouping
        {
        public:
            int idOf(const std::string & name)
            {
                const auto [pos, inserted] =
                  m_Ids.try_emplace(name, static_cast<int>(m_Ids.size()));
                return pos->second;
            }

        private:
            std::map<std::string, int> m_Ids;
        };

        //! Translate one solved element into the estimator's plain-data element, recording its
        //! node coordinates into the shared node table along the way.
        Element toEstimatorElement(const IElementLinear2D & element,
                                   MaterialGrouping & grouping,
                                   std::vector<Point> & nodes)
        {
            const auto nodeIds = element.nodeIndexes();
            const auto gaussCoords = element.gaussPointGlobalCoordinates();
            const auto gaussFlux = element.fluxAtGaussPoints();
            const auto conductivity = element.meanConductivity();
            const double inverse = conductivity == 0.0 ? 0.0 : 1.0 / conductivity;

            Element out;
            out.inverseConstitutive = {inverse, 0.0, 0.0, inverse};
            out.subdomain = grouping.idOf(element.getMaterial().name());

            for(std::size_t vtx = 0; vtx < numOfQuadrilateralNodes; ++vtx)
            {
                const auto nodeId = nodeIds[vtx];
                out.vertexIds[vtx] = static_cast<int>(nodeId);
                if(nodeId >= nodes.size())
                {
                    nodes.resize(nodeId + 1);
                }
                const auto & node = element.getNode(vtx);
                nodes[nodeId] = {node.X(), node.Y()};
            }

            for(std::size_t gpt = 0; gpt < numOfQuadrilateralNodes; ++gpt)
            {
                out.gaussPoints[gpt] = {gaussCoords[gpt][0], gaussCoords[gpt][1]};
                out.flux[gpt] = {gaussFlux[gpt].x, gaussFlux[gpt].y};
            }

            return out;
        }
    }

    DomainEstimate estimateError(const IDomain & domain, const double targetPercent)
    {
        Input input;
        input.targetPercent = targetPercent;

        MaterialGrouping grouping;
        for(const auto & element : domain.elements().elements())
        {
            input.elements.push_back(toEstimatorElement(*element, grouping, input.nodes));
        }

        auto estimated = estimate(input);
        Result result = estimated.has_value() ? std::move(estimated).value() : Result{};
        return DomainEstimate{std::move(input), std::move(result)};
    }
}
