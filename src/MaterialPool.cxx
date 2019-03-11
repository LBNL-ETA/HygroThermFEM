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
    }

    const IMaterial & MaterialPool::createSolidMaterial(
      const std::string & Name,
      double Density,
      double Porosity,
      double HeatCapacity,
      double DiffusionResistanceFactor,
      const std::vector<std::pair<double, double>> & ThermalConductivity,
      const std::vector<std::pair<double, double>> & LiquidTransportCurve,
      const std::vector<std::pair<double, double>> & SorptionCurve,
      double emissivity)
    {
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

    const IMaterial & MaterialPool::createGas(const std::string & Name)
    {
        m_Materials.emplace(std::make_pair(Name, std::unique_ptr<Gas>(new Gas(Name))));
        return *m_Materials.at(Name);
    }

    const IMaterial & MaterialPool::material(const std::string & name) const
    {
        return *m_Materials.at(name);
    }

    std::vector<std::string> MaterialPool::getMaterials(MaterialType materialType) const
    {
        std::vector<std::string> result;
        for(auto & mat : m_Materials)
        {
            if(mat.second->materialType() == materialType)
            {
                result.push_back(mat.second->name());
            }
        }
        return result;
    }

}   // namespace HygroThermFEM
