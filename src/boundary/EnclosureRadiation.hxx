#pragma once

#include <cstddef>
#include <map>
#include <vector>

namespace HygroThermFEM
{
    class Nodes;

    //! \brief How an enclosure surface's temperature enters the radiation exchange.
    enum class EnclosureSurfaceTemperature
    {
        //! Legacy-Conrad convention (bcrad2): each segment is isothermal at the fourth-power mean
        //! of its node temperatures, and the segment's uniform net flux is shared equally by its
        //! two nodes. Reproduces Conrad's enclosure results.
        SegmentIsothermal,
        //! The surface temperature varies along the segment (Gauss-point interpolation in the
        //! boundary condition, linear-mean segment temperature in the radiosity). More refined
        //! than Conrad; produces slightly different corner temperatures on coarse meshes.
        LocalTemperature
    };

    //! \brief What kind of surface an enclosure segment is.
    enum class EnclosureSurfaceKind
    {
        //! On the FEM mesh: coordinates and temperature come from node1/node2 in the node pool, and
        //! the segment gets a finite-element boundary condition.
        Meshed,
        //! Not on the FEM mesh: an emitter at a known temperature with explicit coordinates - e.g.
        //! the IGU gap end-faces in Condensation Index mode, which close the gap cavity at the
        //! computed hole (cavity) temperature. It participates in the radiosity but has no FEM BC.
        FixedTemperature
    };

    //! \brief One radiation surface segment of an enclosure.
    //!
    //! A segment is a two-node boundary line with a surface emissivity, belonging to an enclosure
    //! identified by enclosureId. Segments with the same enclosureId exchange radiation. The kind
    //! selects which fields are used: Meshed uses node1/node2 (node pool); FixedTemperature uses
    //! fixedTemperature and the start/end coordinates.
    struct EnclosureRadiationSegment
    {
        EnclosureSurfaceKind kind{EnclosureSurfaceKind::Meshed};

        std::size_t node1{0};
        std::size_t node2{0};
        double emissivity{0.9};
        std::size_t enclosureId{0};

        double fixedTemperature{0.0};
        double startX{0.0};
        double startY{0.0};
        double endX{0.0};
        double endY{0.0};
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
        //! \param surfaceTemperature Segment surface-temperature convention (Conrad-compatible
        //! isothermal by default).
        EnclosureRadiation(Nodes & nodePool,
                           std::vector<EnclosureRadiationSegment> segments,
                           std::map<std::size_t, double> openEnclosureTemperatures,
                           bool smoothViewFactors = true,
                           EnclosureSurfaceTemperature surfaceTemperature =
                             EnclosureSurfaceTemperature::SegmentIsothermal);

        [[nodiscard]] std::size_t numberOfSegments() const;

        //! Effective radiant temperature [C] seen by the given segment, from the enclosure
        //! radiosity at the current nodal temperatures. Re-solves only when temperatures change.
        [[nodiscard]] double effectiveRadiantTemperature(std::size_t segmentIndex);

        //! The surface-temperature convention this coordinator (and its BCs) follow.
        [[nodiscard]] EnclosureSurfaceTemperature surfaceTemperatureModel() const;

        //! Segment surface temperature [C] under the current convention (for SegmentIsothermal,
        //! the fourth-power mean of the node temperatures - the value the segment radiates at).
        [[nodiscard]] double segmentSurfaceTemperature(std::size_t segmentIndex);

    private:
        void solveIfNeeded();
        [[nodiscard]] std::vector<double> currentSegmentTemperatures() const;

        Nodes & m_Nodes;
        std::vector<EnclosureRadiationSegment> m_Segments;
        EnclosureSurfaceTemperature m_SurfaceTemperature;

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
