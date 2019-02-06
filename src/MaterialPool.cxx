#include <algorithm>

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

    const Material & MaterialPool::createMaterial(
      const std::string & Name,
      const double Density,
      const double Porosity,
      const double HeatCapacity,
      const double DiffusionResistanceFactor,
	  const std::vector<std::pair<double, double>> & ThermalConductivity,
      const std::vector<std::pair<double, double>> & LiquidTransportCurve,
      const std::vector<std::pair<double, double>> & SorptionCurve)
    {
        m_Materials.emplace(std::make_pair(Name,
                                           Material(Name,
                                                    Density,
                                                    Porosity,
                                                    HeatCapacity,
                                                    DiffusionResistanceFactor,
                                                    ThermalConductivity,
                                                    LiquidTransportCurve,
                                                    SorptionCurve)));
        return m_Materials.at(Name);
    }

    const Material & MaterialPool::material(const std::string & name) const
    {
        return m_Materials.at(name);
    }

}   // namespace HygroThermFEM
