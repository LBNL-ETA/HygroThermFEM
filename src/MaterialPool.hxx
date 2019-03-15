#pragma once

#include <map>

#include "Material.hxx"

namespace HygroThermFEM
{
    class MaterialPool
    {
    public:
        static MaterialPool & Instance();

        const IMaterial & createSolidMaterial(const std::string & Name,
                                         double Density,
                                         double Porosity,
                                         double HeatCapacity,
                                         double DiffusionResistanceFactor,
                                         const std::vector<std::pair<double, double>> &
                                         ThermalConductivity,
                                         const std::vector<std::pair<double, double>> &
                                         LiquidTransportCurve,
                                         const std::vector<std::pair<double, double>> &
                                         SorptionCurve,
                                         double emissivity = 0.9);

        const IMaterial & createGas(const std::string & Name);

        const IMaterial & material(const std::string & name) const;

        void clear();

        //! Returns all material names of given type
        std::vector<std::string> getMaterials(MaterialType materialType) const;

    private:
        MaterialPool() = default;

        ~MaterialPool() = default;

        std::map<std::string, std::unique_ptr<IMaterial>> m_Materials;
    };

}   // namespace HygroThermFEM
