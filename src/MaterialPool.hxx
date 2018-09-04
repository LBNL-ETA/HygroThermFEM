#pragma once

#include <map>

#include "Material.hxx"

namespace MoisThermFEM
{
    class MaterialPool
    {
    public:
        static MaterialPool & Instance();

        const Material &
          createMaterial(const std::string & Name,
                         double Density,
                         double Porosity,
                         double HeatCapacity,
                         double ThermalConductivity,
                         double DiffusionResistanceFactor,
                         const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                         const std::vector<std::pair<double, double>> & SorptionCurve);

        const Material & material(const std::string & name) const;

        void clear();

    private:
        MaterialPool() = default;

        ~MaterialPool() = default;

        std::map<std::string, Material> m_Materials;
    };

}   // namespace MoisThermFEM
