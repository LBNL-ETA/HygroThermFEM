#include "MaterialMissingProperties.hxx"

namespace HygroThermFEM
{
    //////////////////////////////////////////////////////////////////////////
    ///  MaterialMissingProperties
    //////////////////////////////////////////////////////////////////////////
    std::vector<std::string> MaterialMissingProperties::missingPropertiesMessage() const
    {
        std::vector<std::string> message{};
        if(isMissingAnyProperty())
        {
            message.push_back("Material " + materialName + " is missing following properties:");
            if(Density)
            {
                message.emplace_back("- Density");
            }
            if(Emissivity)
            {
                message.emplace_back("- Emissivity");
            }
            if(Porosity)
            {
                message.emplace_back("- Porosity");
            }
            if(SpecificHeatCapacityDry)
            {
                message.emplace_back("- Specific Heat Capacity Dry");
            }
            if(ThermalConductivityDry)
            {
                message.emplace_back("- Thermal Conductivity Dry");
            }
            if(WaterVaporDiffusionResistanceFactor)
            {
                message.emplace_back("- Water Vapor Diffusion Resistance Factor");
            }
            if(MoistureStorageFunction)
            {
                message.emplace_back("- Moisture Storage Function");
            }
            if(LiquidTransportationSuction)
            {
                message.emplace_back("- Liquid Transportation Suction Curve");
            }
            if(LiquidTransportationRedistribution)
            {
                message.emplace_back("- Liquid Transportation Redistribution Curve");
            }
            if(ThermalConductivityMoistureAndTemperatureDependent)
            {
                message.emplace_back("- Thermal Conductivity Moisture and Temperature Dependent");
            }
        }
        return message;
    }

    bool MaterialMissingProperties::isMissingAnyProperty() const
    {
        return Density || Emissivity || Porosity || SpecificHeatCapacityDry
               || WaterVaporDiffusionResistanceFactor || MoistureStorageFunction
               || LiquidTransportationSuction || LiquidTransportationRedistribution
               || ThermalConductivityMoistureAndTemperatureDependent;
    }

    bool isMaterialLibraryCorrect(
      const std::vector<MaterialMissingProperties> & missingProperties)
    {
        bool isCorrect{true};
        for(const auto & mat : missingProperties)
        {
            isCorrect = isCorrect && !mat.isMissingAnyProperty();
        }
        return isCorrect;
    }

}   // namespace HygroThermFEM