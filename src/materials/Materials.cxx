#include <stdexcept>
#include <variant>

#include "Materials.hxx"

namespace HygroThermFEM
{
    void Materials::clear()
    {
        m_Pool.clear();
    }

    const IMaterial & Materials::createSolidMaterial(SolidMaterialParams params)
    {
        checkIfMaterialExists(params.name);
        const auto name = params.name;
        m_Pool.try_emplace(name, std::in_place_type<SolidMaterial>, std::move(params));
        return std::get<SolidMaterial>(m_Pool.at(name));
    }

    IMaterial & Materials::createSolidMaterial(std::string Name)
    {
        checkIfMaterialExists(Name);
        const auto name = Name;
        m_Pool.try_emplace(name, std::in_place_type<SolidMaterial>, std::move(Name));
        return std::get<SolidMaterial>(m_Pool.at(name));
    }

    const IGas & Materials::createGas(const std::string & name,
                                      const CavityStandard cavityStandard, Gases::CGas gas)
    {
        // Preserve the historical overwrite-on-recreate behaviour. The variant alternatives carry
        // const members (so they are not move-assignable); erase then construct in place rather
        // than assigning.
        m_Pool.erase(name);
        m_Pool.try_emplace(name, std::in_place_type<Gas>, name, cavityStandard, gas);
        return std::get<Gas>(m_Pool.at(name));
    }

    IMaterial & Materials::material(const std::string & name) const
    {
        return std::visit([](auto & entry) -> IMaterial & { return entry; }, m_Pool.at(name));
    }

    IGas & Materials::gas(const std::string & name) const
    {
        return std::get<Gas>(m_Pool.at(name));
    }

    std::vector<std::string> Materials::getMaterials() const
    {
        std::vector result{getSolidMaterials()};
        const std::vector gases{getGases()};
        result.insert(result.end(), gases.begin(), gases.end());
        return result;
    }

    std::vector<std::string> Materials::getSolidMaterials() const
    {
        std::vector<std::string> result;
        for(const auto & [name, entry] : m_Pool)
        {
            if(std::holds_alternative<SolidMaterial>(entry))
            {
                result.push_back(name);
            }
        }
        return result;
    }

    std::vector<std::string> Materials::getGases() const
    {
        std::vector<std::string> result;
        for(const auto & [name, entry] : m_Pool)
        {
            if(std::holds_alternative<Gas>(entry))
            {
                result.push_back(name);
            }
        }
        return result;
    }

    void Materials::checkIfMaterialExists(const std::string & materialName) const
    {
        const auto found = m_Pool.find(materialName);
        if(found == m_Pool.end())
        {
            return;
        }
        if(std::holds_alternative<Gas>(found->second))
        {
            throw std::runtime_error("Gas with given name is already inserted in model.");
        }
        throw std::runtime_error("Material with given name is already inserted in model.");
    }

}   // namespace HygroThermFEM
