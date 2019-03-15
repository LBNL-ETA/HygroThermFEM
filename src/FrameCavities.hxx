#pragma once

#include "Elements2D.hxx"
#include "KeffCavity.hxx"

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////
    ///  EquivalentFrameCavity
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Used to create equivalent frame cavity.
    //!
    //! Class first performs rectangularization and then applies one of two standards for
    //! calculating equivalent thermal conductivity.
    class EquivalentFrameCavity
    {
    public:
        explicit EquivalentFrameCavity(const std::vector<size_t> & nodes);

        double L() const;
        double H() const;

        double area() const;

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

            //! Emissivity of the segment
            double emissivity() const;

            //! Segment length
            double length() const;

            //! Segment's average temperature
            double averageTemperature() const;

            const Node2D & firstNode() const;

            //! Helper function that is used in area calculation for entire frame cavity.
            double crossCalc() const;

            //! Returns side (Left, Right, Top or Bottom) to which segment belongs.
            const Side & side() const;

        private:
            Side calcSide(const Node2D & node1, const Node2D & node2) const;

            const Node2D & node1;
            const Node2D & node2;
            double m_Emissivity;
            // Length is calculated once in segment constructor. It is just faster to execute.
            const double m_Length;
            const Side m_Side;
        };

        std::string findCommonMaterial(const Node2D & node1, const Node2D & node2) const;

        std::vector<Segment> buildSegments(const std::vector<size_t> & nodes);

        EquivalentFrameCavity::Size calcSize(double area) const;

        const std::vector<Segment> m_Segments;
        const double m_Area;
        const Size m_Size;
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///  EquivalentFrameCavities
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Class to calculate equivalent frame cavities from mesh model
    class EquivalentFrameCavities
    {
    public:
        //! Construction of frame cavities.
        //!
        //! \param elements: All elements from the domain.
        explicit EquivalentFrameCavities(const ElementsLinear2D & elements);

        //! Returns rectrangular frame cavity
        //!
        //! \param frameCavityName Frame cavity name for which geometry will be returned.
        //! \return Equivalent rectangular frame cavity
        EquivalentFrameCavity getCavity(const std::string & frameCavityName) const;

        //! \brief Calculation of eqivalent thermal conductivity.
        //!
        //! \param frameCavityName Frame cavity name for which thermal conductivity will be
        //! calculated. \return Equivalent thermal conductivity value.
        double thermalConductivity(const std::string & frameCavityName) const;

    private:
        //! Function to return boundary nodes for given frame cavity
        //!
        //! \param frameCavityName Frame cavity name for which boundary nodes will be returned
        //! \return Nodes that form frame cavity boundary
        const std::vector<size_t> & boundaryNodes(const std::string & frameCavityName) const;

        //! Helper class used in algorithm to determine frame cavity boundaries
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

        //! Function to perform calculation of equivalent frame cavities over the entire domain.
        void calculateEquivalentFrameCavities();

        //! Calculate edges of frame cavity with lines in no specific order
        //!
        //! \param elNodes All elements with nodes that form single frame cavity.
        //! \return Set of lines that form edge of frame cavity. Lines are not sorted
        static std::set<line> getEdges(const std::vector<std::vector<size_t>> & elNodes);

        //! Calculate nodes in ordered way. These nodes are boundary of the frame cavity. One frame
        //! cavity is processed at the time.
        //!
        //! \param allEdges Set of unordered lines that form frame cavity
        //! \return Sorted vector of nodes that form frame cavity boundary
        static std::vector<size_t> edgeNodesOrdered(std::set<line> & allEdges);

        const ElementsLinear2D & m_Elements;

        //! Keeps boundary nodes for every frame cavity in the domain.
        std::map<std::string, std::vector<size_t>> m_BoundaryNodes;
    };

}   // namespace HygroThermFEM
