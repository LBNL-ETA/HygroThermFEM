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
    ///  EquivalentFrameCavities
    ///////////////////////////////////////////////////////////////////////////////

    //! \brief Class to calculate equivalent frame cavities from mesh model
    class EquivalentFrameCavities
    {
    public:
        explicit EquivalentFrameCavities(const ElementsLinear2D & m_Elements);

    private:
        //! Helper function used in algorithm to determine frame cavity boundaries
        class line
        {
        public:
            line(size_t n1, size_t n2) : n1(n1), n2(n2)
            {}

            size_t getN1() const
            {
                return n1;
            }
            size_t getN2() const
            {
                return n2;
            }

            bool operator<(const line & rhs) const
            {
                if(n1 < rhs.n1)
                    return true;
                if(rhs.n1 < n1)
                    return false;
                return n2 < rhs.n2;
            }

            bool operator>(const line & rhs) const
            {
                return rhs < *this;
            }

            bool operator<=(const line & rhs) const
            {
                return !(rhs < *this);
            }

            bool operator>=(const line & rhs) const
            {
                return !(*this < rhs);
            }

        private:
            size_t n1;
            size_t n2;
        };

        //! Function to perform calculation of equivalent frame cavities over entire domain.
        void calculateEquivalentFrameCavities();

        //! calculate edges of frame cavity with no sorted lines
        static std::set<line> getEdges(const std::vector<std::vector<size_t>> & elNodes);

        //! calculate nodes in ordered way. These nodes are boundary of the frame cavity. One frame
        //! cavity is processed at the time.
        static std::vector<size_t> edgeNodesOrdered(std::set<line> & allEdges);

        std::vector<FrameCavity> m_EquivalentFrameCavities;
        const ElementsLinear2D & m_Elements;
    };

}   // namespace HygroThermFEM
