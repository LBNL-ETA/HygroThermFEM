#include "EnclosureRadiation.hxx"

#include <cmath>
#include <set>

#include <WCEViewer.hpp>

#include "Common.hxx"
#include "Nodes.hxx"


namespace HygroThermFEM
{
    namespace
    {
        constexpr double zeroCelsiusInKelvin = 273.15;
    }

    EnclosureRadiation::EnclosureRadiation(Nodes & nodePool,
                                           std::vector<EnclosureRadiationSegment> segments,
                                           std::map<std::size_t, double> openEnclosureTemperatures,
                                           const bool smoothViewFactors) :
        m_Nodes(nodePool),
        m_Segments(std::move(segments))
    {
        const auto count = m_Segments.size();

        std::set<std::size_t> openIds;
        for(const auto & [enclosureId, temperature] : openEnclosureTemperatures)
        {
            openIds.insert(enclosureId);
        }

        // View factors from the WCE engine, built from the segment node coordinates (geometry is
        // fixed for the whole solve).
        std::vector<Viewer::RadiationSegment> radiationSegments;
        radiationSegments.reserve(count);
        for(const auto & segment : m_Segments)
        {
            if(segment.kind == EnclosureSurfaceKind::FixedTemperature)
            {
                // Non-meshed enclosure surface (e.g. an IGU gap end-face): use the supplied
                // coordinates; its nodes are not in the FEM pool.
                radiationSegments.push_back({.startPoint = {segment.startX, segment.startY},
                                             .endPoint = {segment.endX, segment.endY},
                                             .emissivity = segment.emissivity,
                                             .enclosureId = segment.enclosureId});
            }
            else
            {
                const auto & node1 = m_Nodes.getNode(segment.node1);
                const auto & node2 = m_Nodes.getNode(segment.node2);
                radiationSegments.push_back({.startPoint = {node1.X(), node1.Y()},
                                             .endPoint = {node2.X(), node2.Y()},
                                             .emissivity = segment.emissivity,
                                             .enclosureId = segment.enclosureId});
            }
        }

        const Viewer::ViewFactorOptions options{.leastSquaresSmoothing = smoothViewFactors,
                                                .openEnclosureIds = openIds};
        const auto viewFactors = Viewer::computeEnclosureViewFactors(radiationSegments, {}, options);
        m_ViewFactors = viewFactors.viewFactors.getMatrix();
        m_EnvironmentViewFactor = viewFactors.environmentViewFactor;

        m_Emissivity.reserve(count);
        m_IsOpen.reserve(count);
        m_EnvironmentTempKelvin.assign(count, 0.0);
        for(std::size_t idx = 0; idx < count; ++idx)
        {
            m_Emissivity.push_back(m_Segments[idx].emissivity);
            const bool isOpen = openIds.contains(m_Segments[idx].enclosureId);
            m_IsOpen.push_back(isOpen);
            if(isOpen)
            {
                m_EnvironmentTempKelvin[idx] =
                  celsiusToKelvin(openEnclosureTemperatures.at(m_Segments[idx].enclosureId));
            }
        }
    }

    std::size_t EnclosureRadiation::numberOfSegments() const
    {
        return m_Segments.size();
    }

    std::vector<double> EnclosureRadiation::currentSegmentTemperatures() const
    {
        std::vector<double> result;
        result.reserve(m_Segments.size());
        for(const auto & segment : m_Segments)
        {
            if(segment.kind == EnclosureSurfaceKind::FixedTemperature)
            {
                // Non-meshed surface held at its (computed) fixed temperature, e.g. the IGU gap
                // end-faces at the hole/cavity temperature.
                result.push_back(segment.fixedTemperature);
                continue;
            }
            const auto temperature1 =
              m_Nodes.getNode(segment.node1).property(Variable::temperature);
            const auto temperature2 =
              m_Nodes.getNode(segment.node2).property(Variable::temperature);
            result.push_back(0.5 * (temperature1 + temperature2));
        }
        return result;
    }

    void EnclosureRadiation::solveIfNeeded()
    {
        const auto temperatures = currentSegmentTemperatures();
        if(m_Solved && temperatures == m_CachedTemperatures)
        {
            return;
        }

        const auto count = m_Segments.size();
        const double sigma = Constants::STEFANBOLTZMANN;

        // Blackbody emissive power of each surface and of each open enclosure's environment.
        std::vector<double> emissivePower(count, 0.0);
        std::vector<double> environmentEmissivePower(count, 0.0);
        for(std::size_t idx = 0; idx < count; ++idx)
        {
            const double surfaceKelvin = celsiusToKelvin(temperatures[idx]);
            emissivePower[idx] = sigma * std::pow(surfaceKelvin, 4);
            if(m_IsOpen[idx])
            {
                environmentEmissivePower[idx] = sigma * std::pow(m_EnvironmentTempKelvin[idx], 4);
            }
        }

        // Radiosity system  (I - (1 - eps) F) J = eps Eb + (1 - eps) Fenv Eb_env.
        std::vector<std::vector<double>> systemData(count, std::vector<double>(count, 0.0));
        std::vector<double> rightHandSide(count, 0.0);
        for(std::size_t row = 0; row < count; ++row)
        {
            const double oneMinusEmissivity = 1.0 - m_Emissivity[row];
            for(std::size_t col = 0; col < count; ++col)
            {
                const double identity = (row == col) ? 1.0 : 0.0;
                systemData[row][col] = identity - oneMinusEmissivity * m_ViewFactors[row][col];
            }
            rightHandSide[row] = m_Emissivity[row] * emissivePower[row]
                                 + oneMinusEmissivity * m_EnvironmentViewFactor[row]
                                     * environmentEmissivePower[row];
        }

        const FenestrationCommon::SquareMatrix system{systemData};
        const auto radiosity = system.inverse() * rightHandSide;

        // Irradiation G_i and the effective radiant temperature T_rad,i = (G_i / sigma)^(1/4).
        m_RadiantTemperature.assign(count, 0.0);
        for(std::size_t row = 0; row < count; ++row)
        {
            double irradiation = m_EnvironmentViewFactor[row] * environmentEmissivePower[row];
            for(std::size_t col = 0; col < count; ++col)
            {
                irradiation += m_ViewFactors[row][col] * radiosity[col];
            }
            const double radiantKelvin = std::pow(std::max(irradiation, 0.0) / sigma, 0.25);
            m_RadiantTemperature[row] = radiantKelvin - zeroCelsiusInKelvin;
        }

        m_CachedTemperatures = temperatures;
        m_Solved = true;
    }

    double EnclosureRadiation::effectiveRadiantTemperature(const std::size_t segmentIndex)
    {
        solveIfNeeded();
        return m_RadiantTemperature[segmentIndex];
    }
}   // namespace HygroThermFEM
