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

    //////////////////////////////////////////////////////////////////////////
    ///  MaterialsErrorCheckVector
    //////////////////////////////////////////////////////////////////////////

    MaterialMissingProperties MaterialsErrorCheckVector::operator[](const size_t Index) const
    {
        return m_MaterialMissingProperties[Index];
    }

    bool MaterialsErrorCheckVector::isMaterialLibraryCorrect() const
    {
        bool isCorrect{true};
        for(const auto & mat : m_MaterialMissingProperties)
        {
            isCorrect = isCorrect && !mat.isMissingAnyProperty();
        }
        return isCorrect;
    }

    void MaterialsErrorCheckVector::push_back(MaterialMissingProperties missingProperties)
    {
        m_MaterialMissingProperties.push_back(std::move(missingProperties));
    }

    size_t MaterialsErrorCheckVector::size() const
    {
        return m_MaterialMissingProperties.size();
    }
}   // namespace HygroThermFEM