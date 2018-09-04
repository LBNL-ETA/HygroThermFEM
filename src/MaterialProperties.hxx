#ifndef MOISTHERMFEM_MATERIALPROPERTIES_HXX
#define MOISTHERMFEM_MATERIALPROPERTIES_HXX

#include <memory>

#include "Material.hxx"

namespace MoisThermFEM
{
    class IValue;

    using iValue = std::unique_ptr<IValue>;

    class MaterialProperties
    {
    public:
        static iValue getWaterFill(const Material & mat);

        static iValue getAirFill(const Material & mat);
    };

}   // namespace MoisThermFEM

#endif   // MOISTHERMFEM_MATERIALPROPERTIES_HXX
