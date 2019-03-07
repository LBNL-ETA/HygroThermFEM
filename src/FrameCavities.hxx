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
    //!
    //! Rectangularization is done according to ISO 15099 standard. These calculations are needed
    //! for calculation of thermal conductivity.
    class RectangularizedCavity
    {
    public:
        explicit RectangularizedCavity(const std::vector<size_t> & nodes);

    private:

        //! \brief Local enumerator used to assign segment to certain side of rectangular frame
        //! cavity.
        enum class Side
        {
            Top,
            Bottom,
            Left,
            Right
        };

        struct Size
        {
            double L;
            double H;
        };

        //! \brief Local class segment that keeps basic properties for recangular cavity
        //! calculations.
        class Segment
        {
        public:
            Segment(const Node2D & node1, const Node2D & node2, double emissivity);

            double emissivity() const;
            double length() const;
            double node1Temperature() const;
            double node2Temperature() const;
            double averageTemperature() const;

            const Node2D & firstNode() const;

            //! Helper function that is used in area calculation for entire frame cavity.
            double crossCalc() const;

        private:
            Side calcSide(const Node2D &node1, const Node2D &node2) const;

            const Node2D & node1;
            const Node2D & node2;
            double m_Emissivity;
            // Length is calculated once in segment constructor. It is just faster to execute.
            const double m_Length;
            const Side m_Side;
        };

        std::string findCommonMaterial(const Node2D & node1, const Node2D & node2) const;

        std::vector<Segment> buildSegments(const std::vector<size_t> & nodes);

        double area() const;

        RectangularizedCavity::Size calcSize(double area) const;

        const std::vector<Segment> m_Segments;
        const double m_Area;
        const Size m_Size;
    };

}   // namespace HygroThermFEM
