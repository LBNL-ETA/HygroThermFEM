#pragma once

#include <optional>
#include <functional>

#include <KeffCavity.hxx>
#include "Elements2D.hxx"
#include "EnumerationTemplate.hpp"
#include "Materials.hxx"

namespace HygroThermFEM
{
    class Nodes;
    ///////////////////////////////////////////////////////////////////////////////
    ///  EquivalentGasCavity
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Used to create equivalent properties of air cavity.
    //!
    //! Class first performs rectangularization and then applies one of two standards for
    //! calculating equivalent thermal conductivity.
    class EquivalentGasCavity
    {
    public:
        //! \brief Equivalent frame cavity constructor
        //!
        //! \param nodePool Reference to NodePool for node lookup
        //! \param nodes Nodes from which frame cavity has been created
        //! \param gas Frame cavity is filled with this gas
        //! \param gravityVector Gravity vector in (x, y, z) coordinate system
        explicit EquivalentGasCavity(Nodes & nodePool,
                                     const std::vector<size_t> & nodes,
                                     IGas & gas,
                                     const FenestrationCommon::GravityVector & gravityVector = {
                                         0, -1, 0});

        //! Public function that update gas cavity with new data, calculates new thermal
        void update();

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! \param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

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

        //! \brief Class that will make possible for enumeration to go through for loop in C++ 11 and older standards.
        class EnumSide : public Enum<Side>
        {
        public:
            Iterator begin()
            {
                return Iterator(static_cast<int>(Side::Top));
            }

            Iterator end()
            {
                return Iterator(static_cast<int>(Side::Right) + 1);
            }
        };

        //! \brief Simple structure to hold size of equivalent frame cavity.
        //! L is cavity width
        //! H is cavity height
        struct Size
        {
            double L;
            double H;
        };


        //! \brief Local class segment that keeps basic properties for rectangular cavity
        //! calculations.
        class Segment
        {
        public:
            Segment(const Node2D & node1, const Node2D & node2, double emissivity);

            //! Emissivity of the segment
            [[nodiscard]] double emissivity() const;

            //! Segment length
            [[nodiscard]] double length() const;

            //! Segment's average temperature
            [[nodiscard]] double temperature() const;

            [[nodiscard]] const Node2D & firstNode() const;

            //! Helper function that is used in area calculation for entire frame cavity.
            [[nodiscard]] double crossCalc() const;

            //! Returns side (Left, Right, Top or Bottom) to which segment belongs.
            [[nodiscard]] const Side & side() const;

        private:
            //! Calculates frame cavity side to which segment belongs to.
            //!
            //! \param node1 Segment's first node
            //! \param node2 Segment's second node
            //! \return Side to which segment belongs to.
            [[nodiscard]] Side calcSide(const Node2D & node1, const Node2D & node2) const;

            const Node2D & node1;
            const Node2D & node2;
            double m_Emissivity;
            // Length is calculated once in segment constructor. It is just faster to execute.
            const double m_Length;
            const Side m_Side;
        };

        //! Gets side (temperature and emissivity) based on screen flow direction
        KeffCavity::CavitySide getSide1(KeffCavity::ScreenFlow screenFlow) const;

        //! Gets side (temperature and emissivity) based on screen flow direction
        KeffCavity::CavitySide getSide2(KeffCavity::ScreenFlow screenFlow) const;

        //! Calculates frame cavity area out of segments
        [[nodiscard]] double calcArea() const;

        static std::optional<std::reference_wrapper<const IMaterial>>
            findCommonMaterial(const Node2D & node1, const Node2D & node2);

        static std::vector<Segment> buildSegments(Nodes & nodePool, const std::vector<size_t> & nodes);

        [[nodiscard]] Size calcSize(double area) const;

        //! Group segments at four sides of cavity. This is used in rectangularization algorithm.
        static std::map<Side, std::vector<Segment>>
          groupSegmentSides(const std::vector<Segment> & segments);

        void calcSideEmissivities();

        //! Update temperatures for each side of cavity
        void updateSideTemperatures();

        //! \brief Calculates screen heat flow direction.
        [[nodiscard]] KeffCavity::ScreenFlow heatFlowDirection() const;


        const std::vector<Segment> m_Segments;

        //! Map that keeps segments grouped by sides. This is needed for constant temperatures
        //! recalculation.
        const std::map<Side, std::vector<Segment>> m_SideSegments;

        //! Structure to store every side and its properties.
        std::map<Side, KeffCavity::CavitySide> m_Side;
        const double m_Area;
        const Size m_Size;
        IGas & m_Gas;
        FenestrationCommon::GravityVector m_GravityVector;
    };

    ///////////////////////////////////////////////////////////////////////////////
    ///  EquivalentGasCavities
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Class to calculate equivalent frame cavities from mesh model
    class EquivalentGasCavities
    {
    public:
        //! Construction of frame cavities.
        //!
        //! \param nodePool: Reference to NodePool for node lookup
        //! \param materialPool: Reference to MaterialPool for gas lookups
        //! \param elements: All elements from the domain.
        explicit EquivalentGasCavities(Nodes & nodePool,
                                       Materials & materialPool,
                                       const ElementsLinear2D & elements);

        //! \brief Updates frame cavities with new temperatures.
        void update();

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! \param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

    private:
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
        void createEquivalentFrameCavities();

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

        Nodes & m_NodePool;
        Materials & m_MaterialPool;
        const ElementsLinear2D & m_Elements;

        //! Keeps boundary nodes for every frame cavity in the domain.
        std::map<std::string, std::vector<size_t>> m_BoundaryNodes;
        std::vector<EquivalentGasCavity> m_Cavities;
    };

}   // namespace HygroThermFEM
