#pragma once

#include "Elements2D.hxx"

namespace HygroThermFEM
{
    //! \brief Keeps heat flow direction in the frame cavity relative to the screen
    enum class FrameCavityScreenHeatFlow
    {
        Left,
        Right,
        Up,
        Down
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///  FrameCavity
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Holds additional data calculated in frame cavity
    class FrameCavity
    {
    public:
        FrameCavity(double effectiveConductivity,
                    double horizontalDimension,
                    double verticalDimension,
                    double nusseltNumber,
                    FrameCavityScreenHeatFlow heatFlowDirection,
                    double emissivity1,
                    double temperature1,
                    double emissivity2,
                    double temperature2);

    private:
        struct Side
        {
            double temperature;
            double emissivity;
        };

        double m_EffectiveConductivity;
        double m_HorizontalDimension;
        double m_VerticalDimension;
        double m_NusseltNumber;
        FrameCavityScreenHeatFlow m_HeatFlowDirection;
        Side m_Side1;
        Side m_Side2;
        double m_Emissivity2;
        double m_Temperature2;
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///  FrameCavityBoundaries
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Class to calculate equivalent frame cavities from mesh model
    class FrameCavityBoundaries
    {
    public:
        explicit FrameCavityBoundaries(const ElementsLinear2D & m_Elements);

        const std::vector<size_t> & boundaryNodes(const std::string & frameCavityName) const;

    private:
        //! Helper function used in algorithm to determine frame cavity boundaries
        class line
        {
        public:
            line(size_t n1, size_t n2);

            size_t getN1() const;

            size_t getN2() const;

            bool operator<(const line & rhs) const;

            bool operator>(const line & rhs) const;

            bool operator<=(const line & rhs) const;

            bool operator>=(const line & rhs) const;

        private:
            size_t n1;
            size_t n2;
        };

        //! Function to perform calculation of equivalent frame cavities over entire domain.
        void calculateEquivalentFrameCavities();

        //! Calculate edges of frame cavity with lines in no specific order
        static std::set<line> getEdges(const std::vector<std::vector<size_t>> & elNodes);

        //! Calculate nodes in ordered way. These nodes are boundary of the frame cavity. One frame
        //! cavity is processed at the time.
        static std::vector<size_t> edgeNodesOrdered(std::set<line> & allEdges);

        const ElementsLinear2D & m_Elements;

        //! Keeps boundary nodes for every frame cavity in the domain.
        std::map<std::string, std::vector<size_t>> m_BoundaryNodes;
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///  RectangularizedCavity
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Used to create equivalent rectangular frame cavity.
    class RectangularizedCavity
    {
    public:
        RectangularizedCavity(const std::vector<size_t> & nodes);
    };

}   // namespace HygroThermFEM
