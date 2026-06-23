#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "HygroThermFEM2D.hxx"

namespace TestHelper
{
    //! \brief Description of one horizontal segment of a beam.
    //!
    //! A segment is a single material occupying a contiguous range of element
    //! columns. The segment spans `numElementsX` element columns horizontally
    //! and the entire vertical extent of the beam.
    struct BeamSegment
    {
        std::string material;
        std::size_t numElementsX = 1;
        double width = 0.0;
    };

    //! \brief Builder for creating multi-material beam (wall-like) geometries
    //!        in tests.
    //!
    //! The beam is a 2D rectangular grid of quadrilateral elements with:
    //!   - `numElementsY` element rows in the vertical (y) direction
    //!   - One or more horizontally stacked material segments in the
    //!     x direction
    //!
    //! Node numbering is column-major: starting from the left-most x column at
    //! the bottom and moving up, then advancing to the next x column. This
    //! makes edge node lookup deterministic regardless of segment count.
    //!
    //! Indices used here are 1-based to match the rest of the codebase.
    //!
    //! Example:
    //! \code
    //!   BeamBuilder(multiDomain)
    //!       .xStart(0.0)
    //!       .height(0.05)
    //!       .numElementsY(4)
    //!       .state(state)
    //!       .addSegment({.material = matA.name(), .numElementsX = 5, .width = 0.10})
    //!       .addSegment({.material = matB.name(), .numElementsX = 3, .width = 0.05})
    //!       .build();
    //!
    //!   for(auto [i1, i2] : builder.leftEdge())
    //!       multiDomain.createBC_FixedHc(i1, i2, bcCoeff);
    //! \endcode
    class BeamBuilder
    {
    public:
        explicit BeamBuilder(HygroThermFEM::MultiDomain & domain) : m_Domain(domain)
        {}

        BeamBuilder & xStart(double val)
        {
            m_XStart = val;
            return *this;
        }

        BeamBuilder & height(double val)
        {
            m_Height = val;
            return *this;
        }

        BeamBuilder & numElementsY(std::size_t val)
        {
            m_NumElementsY = val;
            return *this;
        }

        BeamBuilder & state(HygroThermFEM::State val)
        {
            m_State = val;
            return *this;
        }

        BeamBuilder & addSegment(BeamSegment segment)
        {
            m_Segments.push_back(std::move(segment));
            return *this;
        }

        //! Total number of element columns across all segments.
        [[nodiscard]] std::size_t numElementsX() const
        {
            std::size_t total = 0;
            for(const auto & seg : m_Segments)
            {
                total += seg.numElementsX;
            }
            return total;
        }

        //! Number of node columns (one more than element columns).
        [[nodiscard]] std::size_t numNodesX() const
        {
            return numElementsX() + 1;
        }

        //! Number of node rows (one more than element rows).
        [[nodiscard]] std::size_t numNodesY() const
        {
            return m_NumElementsY + 1;
        }

        //! 1-based node index at column `c` (0..numElementsX), row `r` (0..numElementsY).
        //! Column-major: nodes within a column are stored bottom-to-top.
        [[nodiscard]] std::size_t nodeIndex(const std::size_t c, const std::size_t r) const
        {
            return c * numNodesY() + r + 1;
        }

        void build()
        {
            if(m_Segments.empty())
            {
                throw std::runtime_error("BeamBuilder: no segments specified.");
            }
            if(m_NumElementsY == 0)
            {
                throw std::runtime_error("BeamBuilder: numElementsY must be > 0.");
            }

            createNodes();
            createElements();
        }

        //! Returns adjacent node-index pairs along the left edge,
        //! ordered bottom-to-top.
        [[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> leftEdge() const
        {
            std::vector<std::pair<std::size_t, std::size_t>> result;
            result.reserve(m_NumElementsY);
            for(std::size_t r = 0; r < m_NumElementsY; ++r)
            {
                result.emplace_back(nodeIndex(0, r), nodeIndex(0, r + 1));
            }
            return result;
        }

        //! Returns adjacent node-index pairs along the right edge,
        //! ordered bottom-to-top.
        [[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> rightEdge() const
        {
            const std::size_t cRight = numElementsX();
            std::vector<std::pair<std::size_t, std::size_t>> result;
            result.reserve(m_NumElementsY);
            for(std::size_t r = 0; r < m_NumElementsY; ++r)
            {
                result.emplace_back(nodeIndex(cRight, r), nodeIndex(cRight, r + 1));
            }
            return result;
        }

        //! Returns adjacent node-index pairs along the top edge,
        //! ordered left-to-right.
        [[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> topEdge() const
        {
            const std::size_t rTop = m_NumElementsY;
            std::vector<std::pair<std::size_t, std::size_t>> result;
            result.reserve(numElementsX());
            for(std::size_t c = 0; c < numElementsX(); ++c)
            {
                result.emplace_back(nodeIndex(c, rTop), nodeIndex(c + 1, rTop));
            }
            return result;
        }

        //! Returns adjacent node-index pairs along the bottom edge,
        //! ordered left-to-right.
        [[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> bottomEdge() const
        {
            std::vector<std::pair<std::size_t, std::size_t>> result;
            result.reserve(numElementsX());
            for(std::size_t c = 0; c < numElementsX(); ++c)
            {
                result.emplace_back(nodeIndex(c, 0), nodeIndex(c + 1, 0));
            }
            return result;
        }

        //! Edge identifiers for applying boundary conditions.
        enum class Edge { Left, Right, Top, Bottom };

        //! Returns node pairs for the given edge.
        [[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>> edge(const Edge edg) const
        {
            switch(edg)
            {
                case Edge::Left: return leftEdge();
                case Edge::Right: return rightEdge();
                case Edge::Top: return topEdge();
                case Edge::Bottom: return bottomEdge();
            }
            return {};
        }

        // --- Boundary condition helpers ---
        // Each method loops over node pairs on the given edge and forwards
        // to the corresponding MultiDomain::createBC_* overload.

        template<typename Coeff>
        void applyBC_FixedHc(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_FixedHc(idx1, idx2, coeff);
            }
        }

        template<typename Coeff>
        void applyBC_TARPHc(const Edge edg, const Coeff & coeff, const double surfaceTilt = 90)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_TARPHc(idx1, idx2, coeff, surfaceTilt);
            }
        }

        template<typename Coeff>
        void applyBC_ASHRAEInsideHc(const Edge edg,
                                     const Coeff & coeff,
                                     const double surfaceHeight,
                                     const double surfaceTilt = 90)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_ASHRAEInsideHc(idx1, idx2, coeff, surfaceHeight, surfaceTilt);
            }
        }

        template<typename Coeff>
        void applyBC_ASHRAEOutsideHc(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_ASHRAEOutsideHc(idx1, idx2, coeff);
            }
        }

        template<typename Coeff>
        void applyBC_YazdanianKlemsHc(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_YazdanianKlemsHc(idx1, idx2, coeff);
            }
        }

        template<typename Coeff>
        void applyBC_KimuraHc(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_KimuraHc(idx1, idx2, coeff);
            }
        }

        template<typename... Args>
        void applyBC_FixedTemperature(const Edge edg, Args &&... args)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_FixedTemperature(idx1, idx2, std::forward<Args>(args)...);
            }
        }

        template<typename Coeff>
        void applyBC_FixedTemperatureAndHumidity(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_FixedTemperatureAndHumidity(idx1, idx2, coeff);
            }
        }

        template<typename Flux>
        void applyBC_FixedHeatFlux(const Edge edg, Flux && flux)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_FixedHeatFlux(idx1, idx2, std::forward<Flux>(flux));
            }
        }

        void applyBC_BlackBodyRadiation(const Edge edg,
                                         const double emissivity,
                                         const double radiationTemperature)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_BlackBodyRadiation(idx1, idx2, emissivity, radiationTemperature);
            }
        }

        template<typename Coeff>
        void applyBC_BlackBodyRadiation(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_BlackBodyRadiation(idx1, idx2, coeff);
            }
        }

        template<typename Coeff>
        void applyBC_LinearizedRadiation(const Edge edg, const Coeff & coeff)
        {
            for(const auto & [idx1, idx2] : edge(edg))
            {
                m_Domain.createBC_LinearizedRadiation(idx1, idx2, coeff);
            }
        }

    private:
        //! Builds the x-coordinate of every node column by walking segments.
        [[nodiscard]] std::vector<double> buildXCoordinates() const
        {
            std::vector<double> xs;
            xs.reserve(numNodesX());

            double currentX = m_XStart;
            xs.push_back(currentX);
            for(const auto & seg : m_Segments)
            {
                if(seg.numElementsX == 0)
                {
                    throw std::runtime_error("BeamBuilder: segment numElementsX must be > 0.");
                }
                const double dx = seg.width / static_cast<double>(seg.numElementsX);
                for(std::size_t i = 0; i < seg.numElementsX; ++i)
                {
                    currentX += dx;
                    xs.push_back(currentX);
                }
            }
            return xs;
        }

        void createNodes() const
        {
            const auto xs = buildXCoordinates();
            const double dy = m_Height / static_cast<double>(m_NumElementsY);

            // Column-major: for each x column, create nodes bottom-to-top.
            for(std::size_t c = 0; c < xs.size(); ++c)
            {
                for(std::size_t r = 0; r <= m_NumElementsY; ++r)
                {
                    const double y = static_cast<double>(r) * dy;
                    m_Domain.nodes().createNode({.index = nodeIndex(c, r),
                                                 .x = xs[c],
                                                 .y = y,
                                                 .state = m_State});
                }
            }
        }

        void createElements() const
        {
            // Walk segments left-to-right; for each segment, create
            // numElementsX * numElementsY elements with the segment's material.
            std::size_t cBase = 0;
            for(const auto & seg : m_Segments)
            {
                for(std::size_t i = 0; i < seg.numElementsX; ++i)
                {
                    const std::size_t c = cBase + i;
                    for(std::size_t r = 0; r < m_NumElementsY; ++r)
                    {
                        const std::size_t botLeft = nodeIndex(c, r);
                        const std::size_t topLeft = nodeIndex(c, r + 1);
                        const std::size_t botRight = nodeIndex(c + 1, r);
                        const std::size_t topRight = nodeIndex(c + 1, r + 1);

                        // CCW from bottom-left: BL -> BR -> TR -> TL
                        m_Domain.createElement({.node1 = botLeft,
                                                .node2 = botRight,
                                                .node3 = topRight,
                                                .node4 = topLeft,
                                                .material = seg.material});
                    }
                }
                cBase += seg.numElementsX;
            }
        }

        HygroThermFEM::MultiDomain & m_Domain;
        std::vector<BeamSegment> m_Segments;
        double m_XStart = 0.0;
        double m_Height = 0.05;
        std::size_t m_NumElementsY = 1;
        HygroThermFEM::State m_State = {};
    };

}   // namespace TestHelper
