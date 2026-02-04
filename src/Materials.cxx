#include <stdexcept>

#include "Materials.hxx"

namespace HygroThermFEM
{
    void Materials::clear()
    {
        m_Materials.clear();
        m_Gases.clear();
    }

    const IMaterial & Materials::createSolidMaterial(
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
        m_Materials.emplace(Name,
                            std::make_unique<SolidMaterial>(Name,
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
                                                            emissivity));
        return *m_Materials.at(Name);
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
        std::vector<std::string> result{getSolidMaterials()};
        std::vector<std::string> gases{getGases()};
        result.insert(result.end(), gases.begin(), gases.end());
        return result;
    }

    std::vector<std::string> Materials::getSolidMaterials() const
    {
        std::vector<std::string> result;
        for(auto & mat : m_Materials)
        {
            result.push_back(mat.second->name());
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
