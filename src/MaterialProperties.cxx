#include "MaterialProperties.hxx"
#include "State.hxx"

namespace MoisThermFEM
{
    iValue MoisThermFEM::MaterialProperties::getLiquidWaterFill(
		const MoisThermFEM::Material & mat )
    {
        /// Calculate air and water content
        iValue waterContent = TabularFunction::create(mat.sorptionCurve(), Property::humidity);

        /// Calls sorption curve at 100% humidity to get maximum water content
        const auto maxWaterContent = waterContent->value( State( 0, 1, 0, 0 ) );

        return mat.porosity() / maxWaterContent * waterContent;
    }

    iValue MaterialProperties::getAirFill(const Material & mat)
    {
        iValue waterFill = getLiquidWaterFill( mat );
        return mat.porosity() - waterFill;
    }

}   // namespace MoisThermFEM