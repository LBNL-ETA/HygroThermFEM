#include "MaterialPool.hxx"

namespace HygroThermFEM
{
    MaterialPool & HygroThermFEM::MaterialPool::Instance()
    {
        static MaterialPool m_Instance;
        return m_Instance;
    }

    void MaterialPool::clear()
    {
        m_Materials.clear();
        m_Gases.clear();
    }

    const IMaterial & MaterialPool::createSolidMaterial(
            const std::string & Name,
            double Density,
            double Porosity,
            double HeatCapacity,
            double DiffusionResistanceFactor,
            const std::vector<FenestrationCommon::point> &ThermalConductivity,
            const std::vector<FenestrationCommon::point> &LiquidTransportCurve,
            const std::vector<FenestrationCommon::point> &SorptionCurve,
            double emissivity)
    {
        checkIfMaterialExists(Name);
        m_Materials.emplace(
          std::make_pair(Name,
                         std::unique_ptr<SolidMaterial>(new SolidMaterial(Name,
                                                                          Density,
                                                                          Porosity,
                                                                          HeatCapacity,
                                                                          DiffusionResistanceFactor,
                                                                          ThermalConductivity,
                                                                          LiquidTransportCurve,
                                                                          SorptionCurve,
                                                                          emissivity))));
        return *m_Materials.at(Name);
    }

    const IGas & MaterialPool::createGas(const std::string & name,
                                         const CavityStandard cavityStandard)
    {
        m_Gases[name] = std::make_unique<Gas>(name, cavityStandard);
        return *m_Gases.at(name);
    }

    const IMaterial & MaterialPool::material(const std::string & name) const
    {
        if(m_Gases.find(name) != m_Gases.end())
        {
            return *m_Gases.at(name);
        }
        return *m_Materials.at(name);
    }

    IGas & MaterialPool::gas(const std::string & name) const 
    {
        return *m_Gases.at(name);
    }

    std::vector<std::string> MaterialPool::getMaterials() const
    {
        std::vector<std::string> result{getSolidMaterials()};
        std::vector<std::string> gases{getGases()};
        result.insert(result.end(), gases.begin(), gases.end());
        return result;
    }

    std::vector<std::string> MaterialPool::getSolidMaterials() const
    {
        std::vector<std::string> result;
        for(auto & mat : m_Materials)
        {
            result.push_back(mat.second->name());
        }
        return result;
    }

    std::vector<std::string> MaterialPool::getGases() const
    {
        std::vector<std::string> result;
        for(auto & gas : m_Gases)
        {
            result.push_back(gas.second->name());
        }
        return result;
    }

    void MaterialPool::checkIfMaterialExists(const std::string & materialName) const
    {
        if(m_Gases.find(materialName) != m_Gases.end())
        {
            throw std::runtime_error("Gas with given name is already inserted in model.");
        }
        if(m_Materials.find(materialName) != m_Materials.end())
        {
            throw std::runtime_error("Material with given name is already inserted in model.");
        }
    }

}   // namespace HygroThermFEM
