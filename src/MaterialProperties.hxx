#ifndef MOISTHERMFEM_MATERIALPROPERTIES_HXX
#define MOISTHERMFEM_MATERIALPROPERTIES_HXX

#include <memory>

#include "Material.hxx"
#include "State.hxx"
#include "Functions.hxx"

namespace MoisThermFEM
{
    /// Function to return total water content in any state (ice, liquid)
    inline auto getWaterFill(const Material & mat)
      -> decltype(Constant(0) * TabularFunction(mat.sorptionCurve(), Property::humidity))
    {
        /// Calculate air and water content
        auto waterContent = TabularFunction(mat.sorptionCurve(), Property::humidity);

        /// Calls sorption curve at 100% humidity to get maximum water content
        const auto maxWaterContent = waterContent.value(State(0, 1, 0, 0));

        return mat.porosity() / maxWaterContent * waterContent;
    }

    inline auto getLiquidWaterFill(const Material & mat)
      -> decltype(StateValue(Property::liquidPercent) * getWaterFill(mat))
    {
        return StateValue(Property::liquidPercent) * getWaterFill(mat);
    }

    inline auto getIceFill(const Material & mat)
      -> decltype((Constant(1) - StateValue(Property::liquidPercent)) * getWaterFill(mat))
    {
        return (Constant(1) - StateValue(Property::liquidPercent)) * getWaterFill(mat);
    }

    inline auto getAirFill(const Material & mat) -> decltype(mat.porosity() - getWaterFill(mat))
    {
        const auto waterFill = getWaterFill(mat);
        return mat.porosity() - waterFill;
    }

}   // namespace MoisThermFEM

#endif   // MOISTHERMFEM_MATERIALPROPERTIES_HXX
