#ifndef MOISTHERMFEM_MATERIALPROPERTIES_HXX
#define MOISTHERMFEM_MATERIALPROPERTIES_HXX

#include <memory>

#include "Material.hxx"
#include "State.hxx"
#include "Functions.hxx"

namespace MoisThermFEM
{
    inline IOperation<Constant, TabularFunction>
        getLiquidWaterFill(const Material & mat)
    {
        /// Calculate air and water content
        auto waterContent = TabularFunction(mat.sorptionCurve(), Property::humidity);

        /// Calls sorption curve at 100% humidity to get maximum water content
        const auto maxWaterContent = waterContent.value(State(0, 1, 0, 0));

        return mat.porosity() / maxWaterContent * waterContent;
    }

    inline IOperation<Constant, IOperation<Constant, TabularFunction>>
        getAirFill(const Material & mat)
    {
        const auto waterFill = getLiquidWaterFill(mat);
        return mat.porosity() - waterFill;
    }

}   // namespace MoisThermFEM

#endif   // MOISTHERMFEM_MATERIALPROPERTIES_HXX
