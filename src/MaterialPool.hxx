#pragma once

#include <map>

#include "Material.hxx"

namespace HygroThermFEM
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
                         double DiffusionResistanceFactor,
                         const std::vector<std::pair<double, double>> & ThermalConductivity,
                         const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                         const std::vector<std::pair<double, double>> & SorptionCurve,
                         double emissivity = 0.9,
                         MaterialType materialType = MaterialType::Solid);

        const Material & material(const std::string & name) const;

        void clear();

    private:
        MaterialPool() = default;

        ~MaterialPool() = default;

        std::map<std::string, Material> m_Materials;
    };

}   // namespace HygroThermFEM
