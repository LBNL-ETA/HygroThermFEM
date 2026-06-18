#pragma once

#include <vector>

#include "MaterialMissingProperties.hxx"

namespace HygroThermFEM
{
    class MultiDomain;
    class IMaterial;

    //////////////////////////////////////////////////////////////////////////
    ///  MaterialDataChecker
    //////////////////////////////////////////////////////////////////////////

    class MaterialDataChecker
    {
    public:
        explicit MaterialDataChecker(const MultiDomain & multiDomain);

        //! \brief Check all material properties against current engine settings.
        //!
        //! \param isTransientSimulation Different checks are needed for different simulation types.
        //! \return All missing properties for every material.
        MaterialsErrorCheckVector checkMaterialProperties(bool isTransientSimulation);

    private:
        //! \brief Check single material properties against current engine settings
        //!
        //! \param material Material that needs to be checked.
        //! \param isTransientSimulation True if simulation is transient, False if steady-state.
        //! \return Missing properties for given material.
        MaterialMissingProperties checkMaterial(const IMaterial & material,
                                                bool isTransientSimulation) const;

        const MultiDomain & multiDomain;
    };
}   // namespace HygroThermFEM