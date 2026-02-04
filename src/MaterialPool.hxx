#pragma once

#include <map>

#include "Material.hxx"

namespace HygroThermFEM
{
    class MaterialPool
    {
    public:

        const IMaterial & createSolidMaterial(
          const std::string & Name,
          double ThermalConductivityDry,
          double Density,
          double Porosity,
          double HeatCapacity,
          double DiffusionResistanceFactor,
          const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
          double moistureDependentMeasurementTemperature,
          const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
          double temperatureDependentMeasurementHumidity,
          const std::vector<FenestrationCommon::point> & LiquidTransportCurve,
          const std::vector<FenestrationCommon::point> & SorptionCurve,
          double emissivity = 0.9);

        IMaterial & createSolidMaterial(std::string Name);

        const IGas & createGas(const std::string & name,
                               CavityStandard cavityStandard = CavityStandard::ISO15099, Gases::CGas gas = Gases::CGas());

        IMaterial & material(const std::string & name) const;

        IGas & gas(const std::string & name) const;

        void clear();

        //! Returns all material names
        [[nodiscard]] std::vector<std::string> getMaterials() const;

        [[nodiscard]] std::vector<std::string> getSolidMaterials() const;

        [[nodiscard]] std::vector<std::string> getGases() const;

        // Public constructor/destructor for ownership by MultiDomain
        MaterialPool() = default;
        ~MaterialPool() = default;

        // Non-copyable, non-movable
        MaterialPool(const MaterialPool &) = delete;
        MaterialPool & operator=(const MaterialPool &) = delete;
        MaterialPool(MaterialPool &&) = delete;
        MaterialPool & operator=(MaterialPool &&) = delete;

    private:

        //! Important check if material has already been created. Program should not allow user to
        //! create two materials with same name because first material will be overwritten and
        //! incorrect results would be calculated by the program.
        void checkIfMaterialExists(const std::string & materialName) const;

        std::map<std::string, std::unique_ptr<IMaterial>> m_Materials;
        std::map<std::string, std::unique_ptr<IGas>> m_Gases;
    };

}   // namespace HygroThermFEM
