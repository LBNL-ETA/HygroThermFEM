
#include "FrameCavities.hxx"

namespace HygroThermFEM
{
    FrameCavity::FrameCavity(const double effectiveConductivity,
                             const double horizontalDimension,
                             const double verticalDimension,
                             const double nusseltNumber,
                             const HygroThermFEM::FrameCavityScreenHeatFlow heatFlowDirection,
                             const double emissivity1,
                             const double temperature1,
                             const double emissivity2,
                             const double temperature2) :
        m_EffectiveConductivity(effectiveConductivity),
        m_HorizontalDimension(horizontalDimension),
        m_VerticalDimension(verticalDimension),
        m_NusseltNumber(nusseltNumber),
        m_HeatFlowDirection(heatFlowDirection),
        m_Side1{temperature1, emissivity1},
        m_Side2{temperature2, emissivity2}
    {}

    EquivalentFrameCavities::EquivalentFrameCavities(const ElementsLinear2D & m_Elements) :
        m_Elements(m_Elements)
    {}
}   // namespace HygroThermFEM
