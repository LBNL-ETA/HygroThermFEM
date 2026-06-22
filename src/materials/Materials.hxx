#pragma once

#include <map>
#include <variant>

#include "Material.hxx"

namespace HygroThermFEM
{
    class Materials
    {
    public:
        //! \brief Create SolidMaterial using params struct (C++20 designated initializers)
        const IMaterial & createSolidMaterial(SolidMaterialParams params);

        //! \brief Create SolidMaterial with just a name (properties set later via setters)
        IMaterial & createSolidMaterial(std::string Name);

        const IGas & createGas(const std::string & name,
                               CavityStandard cavityStandard = CavityStandard::ISO15099, Gases::CGas gas = Gases::CGas());

        IMaterial & material(const std::string & name) const;

        IGas & gas(const std::string & name) const;

        void clear();

        //! Returns all material names
        [[nodiscard]] std::vector<std::string> getMaterials() const;

        [[nodiscard]] std::vector<std::string> getSolidMaterials() const;

        [[nodiscard]] std::vector<std::string> getGases() const;

        // Public constructor/destructor for ownership by MultiDomain
        Materials() = default;
        ~Materials() = default;

        // Non-copyable, non-movable
        Materials(const Materials &) = delete;
        Materials & operator=(const Materials &) = delete;
        Materials(Materials &&) = delete;
        Materials & operator=(Materials &&) = delete;

    private:

        //! Important check if material has already been created. Program should not allow user to
        //! create two materials with same name because first material will be overwritten and
        //! incorrect results would be calculated by the program.
        void checkIfMaterialExists(const std::string & materialName) const;

        //! A single pool holds both solid materials and gases by value, keyed by name. This
        //! replaces the historical parallel m_Materials / m_Gases maps; gas-only lookups use
        //! std::get<Gas> / std::holds_alternative instead of a separate container. It is mutable
        //! because materials are updated in place during the solve (e.g. a gas's thermal
        //! conductivity), which the engine does even through a const Materials reference.
        using MaterialEntry = std::variant<SolidMaterial, Gas>;
        mutable std::map<std::string, MaterialEntry> m_Pool;
    };

}   // namespace HygroThermFEM
