#include <stdexcept>

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
      const double ThermalConductivityDry,
      const double Density,
      const double Porosity,
      const double HeatCapacity,
      const double DiffusionResistanceFactor,
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      const double moistureDependentMeasurementTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      const double temperatureDependentMeasurementHumidity,
      const std::vector<FenestrationCommon::point> & LiquidTransportCurve,
      const std::vector<FenestrationCommon::point> & SorptionCurve,
      const double emissivity)
    {
        checkIfMaterialExists(Name);
        m_Materials.emplace(std::make_pair(
          Name,
          std::unique_ptr<SolidMaterial>(new SolidMaterial(Name,
                                                           ThermalConductivityDry,
                                                           Density,
                                                           Porosity,
                                                           HeatCapacity,
                                                           DiffusionResistanceFactor,
                                                           thermalConductivityMoistureDependent,
                                                           moistureDependentMeasurementTemperature,
                                                           thermalConductivityTemperatureDependent,
                                                           temperatureDependentMeasurementHumidity,
                                                           LiquidTransportCurve,
                                                           SorptionCurve,
                                                           emissivity))));
        return *m_Materials.at(Name);
    }

    IMaterial & MaterialPool::createSolidMaterial(
      std::string Name)
    {
        checkIfMaterialExists(Name);
        m_Materials.emplace(std::make_pair(Name,new SolidMaterial(Name)));
        return *m_Materials.at(Name);
    }

    const IGas & MaterialPool::createGas(const std::string & name,
                                         const CavityStandard cavityStandard, Gases::CGas gas)
    {
        m_Gases[name] = std::make_unique<Gas>(name, cavityStandard, gas);
        return *m_Gases.at(name);
    }

    IMaterial & MaterialPool::material(const std::string & name) const
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
