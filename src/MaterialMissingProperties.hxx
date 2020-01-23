#pragma once

#include <string>
#include <vector>

namespace HygroThermFEM
{
    //! \brief Hold data for missing properties in single material. This is part of error check.
    struct MaterialMissingProperties
    {
        std::string materialName;
        bool Density{false};
        bool Emissivity{false};
        bool Porosity{false};
        bool SpecificHeatCapacityDry{false};
        bool ThermalConductivityDry{false};
        bool WaterVaporDiffusionResistanceFactor{false};
        bool MoistureStorageFunction{false};
        bool LiquidTransportationSuction{false};
        bool LiquidTransportationRedistribution{false};
        bool ThermalConductivityMoistureAndTemperatureDependent{false};

        //! \brief Returns human-readable message about missing properties
        [[nodiscard]] std::vector<std::string> missingPropertiesMessage() const;

        bool isMissingAnyProperty() const;
    };

    //////////////////////////////////////////////////////////////////////////
    ///  MaterialsErrorCheckVector
    //////////////////////////////////////////////////////////////////////////

    //! \brief Simple class that holds info about error data in materials.
    class MaterialsErrorCheckVector
    {
    public:
        MaterialMissingProperties operator[](const size_t Index) const;

        [[nodiscard]] bool isMaterialLibraryCorrect() const;

        void push_back(MaterialMissingProperties missingProperties);

        size_t size() const;

    private:
        std::vector<MaterialMissingProperties> m_MaterialMissingProperties;
    };
}   // namespace HygroThermFEM