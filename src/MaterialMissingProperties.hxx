#pragma once

#include <string>

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
        [[nodiscard]] std::string missingPropertiesMessage() const;

        bool isMissingAnyProperty() const;
    };
}   // namespace HygroThermFEM