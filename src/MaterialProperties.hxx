#ifndef MOISTHERMFEM_MATERIALPROPERTIES_HXX
#define MOISTHERMFEM_MATERIALPROPERTIES_HXX

#include <memory>

#include "Material.hxx"

namespace MoisThermFEM {

	class IValue;

	class MaterialProperties {
	public:
		static std::shared_ptr< MoisThermFEM::IValue > getWaterFill( const Material & mat );

		static std::shared_ptr< MoisThermFEM::IValue > getAirFill( const Material & mat );

	};

}

#endif //MOISTHERMFEM_MATERIALPROPERTIES_HXX
