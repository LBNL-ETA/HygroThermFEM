#include <stdexcept>

#include "Materials.hxx"

#include <ranges>

namespace HygroThermFEM
{
    void Materials::clear()
    {
        m_Materials.clear();
        m_Gases.clear();
    }

    const IMaterial & Materials::createSolidMaterial(SolidMaterialParams params)
    {
        checkIfMaterialExists(params.name);
        const auto name = params.name;
        m_Materials.emplace(name, std::make_unique<SolidMaterial>(std::move(params)));
        return *m_Materials.at(name);
    }

    IMaterial & Materials::createSolidMaterial(std::string Name)
    {
        checkIfMaterialExists(Name);
        const auto name = Name;
        m_Materials.emplace(name, std::make_unique<SolidMaterial>(std::move(Name)));
        return *m_Materials.at(name);
    }

    const IGas & Materials::createGas(const std::string & name,
                                         const CavityStandard cavityStandard, Gases::CGas gas)
    {
        m_Gases[name] = std::make_unique<Gas>(name, cavityStandard, gas);
        return *m_Gases.at(name);
    }

    IMaterial & Materials::material(const std::string & name) const
    {
        if(m_Gases.contains(name))
        {
            return *m_Gases.at(name);
        }
        return *m_Materials.at(name);
    }

    IGas & Materials::gas(const std::string & name) const
    {
        return *m_Gases.at(name);
    }

    std::vector<std::string> Materials::getMaterials() const
    {
        std::vector result{getSolidMaterials()};
        std::vector gases{getGases()};
        result.insert(result.end(), gases.begin(), gases.end());
        return result;
    }

    std::vector<std::string> Materials::getSolidMaterials() const
    {
        std::vector<std::string> result;
        for(const auto & val : m_Materials | std::views::values)
        {
            result.push_back(val->name());
        }
        return result;
    }

    std::vector<std::string> Materials::getGases() const
    {
        std::vector<std::string> result;
        for(auto & gas : m_Gases)
        {
            result.push_back(gas.second->name());
        }
        return result;
    }

    void Materials::checkIfMaterialExists(const std::string & materialName) const
    {
        if(m_Gases.contains(materialName))
        {
            throw std::runtime_error("Gas with given name is already inserted in model.");
        }
        if(m_Materials.contains(materialName))
        {
            throw std::runtime_error("Material with given name is already inserted in model.");
        }
    }

}   // namespace HygroThermFEM
