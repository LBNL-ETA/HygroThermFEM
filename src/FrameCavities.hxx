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

    //! \brief Class to calculate equivalent frame cavities from mesh model
    class EquivalentFrameCavities
    {
    public:
        EquivalentFrameCavities(const ElementsLinear2D &m_Elements);

    private:
        const ElementsLinear2D & m_Elements;
    };

}   // namespace HygroThermFEM