#pragma once

#include <memory>

#include "Material.hxx"
#include "State.hxx"
#include "Functions.hxx"

namespace MoisThermFEM
{
    /// Function to return total water content in liquid + ice
    inline auto getMaterialWaterFill( const Material & mat )
      -> decltype(1.0 * TabularFunction(mat.sorptionCurve(), Property::humidity))
    {
        /// Calculate air and water content
        auto waterContent = TabularFunction(mat.sorptionCurve(), Property::humidity);

        const auto maxWaterContent = waterContent.max();

        return mat.porosity() / maxWaterContent * waterContent;
    }

    /// Water content in liquid state
    inline auto getMaterialLiquidWaterFill( const Material & mat )
      -> decltype(StateValue(Property::liquidPercent) * getMaterialWaterFill( mat ))
    {
        return StateValue(Property::liquidPercent) * getMaterialWaterFill( mat );
    }

    /// Water content in frozen state
    inline auto getMaterialIceFill( const Material & mat )
      -> decltype((1.0 - StateValue(Property::liquidPercent)) * getMaterialWaterFill( mat ))
    {
        return (1.0 - StateValue(Property::liquidPercent)) * getMaterialWaterFill( mat );
    }

    /// Air content
    inline auto getMaterialAirFill( const Material & mat ) -> decltype(mat.porosity() -
		getMaterialWaterFill( mat ))
    {
        const auto waterFill = getMaterialWaterFill( mat );
        return mat.porosity() - waterFill;
    }

}   // namespace MoisThermFEM