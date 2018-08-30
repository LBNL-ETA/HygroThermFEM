#include "MaterialProperties.hxx"
#include "State.hxx"

namespace MoisThermFEM {

	using iValue = std::shared_ptr< MoisThermFEM::IValue >;

	iValue MoisThermFEM::MaterialProperties::getWaterFill( const MoisThermFEM::Material & mat ) {
		/// Calculate air and water content
		iValue waterContent(
				std::make_shared< MoisThermFEM::TabularFunction >( mat.sorptionCurve(),
																													 Property::humidity ) );

		/// Calls sorption curve at 100% humidity to get maximum water content
		auto maxWaterContent = waterContent->value( State( 0, 1, 0 ) );

		return mat.porosity() / maxWaterContent * waterContent;
	}

	iValue MaterialProperties::getAirFill( const Material & mat ) {
		iValue waterFill = getWaterFill( mat );
		return mat.porosity() - waterFill;
	}

}