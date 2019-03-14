#pragma once

#include <map>

#include "Material.hxx"

namespace HygroThermFEM
{
    class MaterialPool
    {
    public:
        static MaterialPool & Instance();

        const IMaterial &
          createSolidMaterial(const std::string & name,
                              double Density,
                              double Porosity,
                              double HeatCapacity,
                              double DiffusionResistanceFactor,
                              const std::vector<std::pair<double, double>> & ThermalConductivity,
                              const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                              const std::vector<std::pair<double, double>> & SorptionCurve,
                              double emissivity = 0.9);

        const IGas & createGas(const std::string & name,
                               const CavityStandard cavityStandard = CavityStandard::ISO15099);

        const IMaterial & material(const std::string & name) const;

        void clear();

        //! Returns all material names
        std::vector<std::string> getMaterials() const;

        std::vector<std::string> getSolidMaterials() const;

        std::vector<std::string> getGases() const;

    private:
        MaterialPool() = default;

        ~MaterialPool() = default;

        //! Important check if material has already been created. Program should not allow user to
        //! create two materials with same name because first material will be overwritten and
        //! incorrect results would be calculated by the program.
        void checkIfMaterialExists(const std::string & materialName) const;

        std::map<std::string, std::unique_ptr<IMaterial>> m_Materials;
        std::map<std::string, std::unique_ptr<IGas>> m_Gases;
    };

}   // namespace HygroThermFEM
