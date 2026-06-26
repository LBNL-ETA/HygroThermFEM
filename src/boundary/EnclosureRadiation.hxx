#pragma once

#include <cstddef>
#include <map>
#include <vector>

namespace HygroThermFEM
{
    class Nodes;

    //! \brief One radiation surface segment of an enclosure.
    //!
    //! A segment is a two-node boundary line with a surface emissivity, belonging to an enclosure
    //! identified by enclosureId. Segments with the same enclosureId exchange radiation.
    struct EnclosureRadiationSegment
    {
        std::size_t node1{0};
        std::size_t node2{0};
        double emissivity{0.9};
        std::size_t enclosureId{0};
    };

    //! \brief Grey-enclosure radiosity solver.
    //!
    //! Gives each enclosure segment the effective radiant temperature of its surroundings, so a
    //! per-segment radiation boundary condition can use the same linearized form as
    //! BlackBodyRadiationBC (just with this dynamic radiant temperature instead of a fixed one).
    //! View factors come from the WCE engine (Viewer::computeEnclosureViewFactors), computed once
    //! at construction since the geometry is fixed; the radiosity solve is memoized on the current
    //! segment temperatures, so the first segment query in an assembly pass solves and the rest of
    //! the pass reuse the result.
    class EnclosureRadiation
    {
    public:
        //! \param nodePool Node pool to read coordinates (at setup) and temperatures (each solve).
        //! \param segments Enclosure radiation segments.
        //! \param openEnclosureTemperatures For each OPEN (auto) enclosure id, the environment
        //! temperature [C]. Enclosures absent from the map are closed.
        //! \param smoothViewFactors Apply least-squares smoothing to closed enclosures.
        EnclosureRadiation(Nodes & nodePool,
                           std::vector<EnclosureRadiationSegment> segments,
                           std::map<std::size_t, double> openEnclosureTemperatures,
                           bool smoothViewFactors = true);

        [[nodiscard]] std::size_t numberOfSegments() const;

        //! Effective radiant temperature [C] seen by the given segment, from the enclosure
        //! radiosity at the current nodal temperatures. Re-solves only when temperatures change.
        [[nodiscard]] double effectiveRadiantTemperature(std::size_t segmentIndex);

    private:
        void solveIfNeeded();
        [[nodiscard]] std::vector<double> currentSegmentTemperatures() const;

        Nodes & m_Nodes;
        std::vector<EnclosureRadiationSegment> m_Segments;

        std::vector<std::vector<double>> m_ViewFactors;   //!< block-diagonal, [i][j]
        std::vector<double> m_EnvironmentViewFactor;      //!< per segment (open-enclosure deficit)
        std::vector<double> m_Emissivity;                 //!< per segment
        std::vector<bool> m_IsOpen;                       //!< per segment: in an open enclosure
        std::vector<double> m_EnvironmentTempKelvin;      //!< per segment (open only)

        std::vector<double> m_RadiantTemperature;         //!< cached effective radiant temp [C]
        std::vector<double> m_CachedTemperatures;         //!< segment temps the cache was solved for
        bool m_Solved{false};
    };
}   // namespace HygroThermFEM
