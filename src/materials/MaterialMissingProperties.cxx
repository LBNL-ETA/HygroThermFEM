#include "MaterialMissingProperties.hxx"

namespace HygroThermFEM
{
    //////////////////////////////////////////////////////////////////////////
    ///  MaterialMissingProperties
    //////////////////////////////////////////////////////////////////////////
    std::vector<std::string> MaterialMissingProperties::missingPropertyNames() const
    {
        std::vector<std::string> names;
        if(Density)
        {
            names.emplace_back("Density");
        }
        if(Emissivity)
        {
            names.emplace_back("Emissivity");
        }
        if(Porosity)
        {
            names.emplace_back("Porosity");
        }
        if(SpecificHeatCapacityDry)
        {
            names.emplace_back("Specific Heat Capacity Dry");
        }
        if(ThermalConductivityDry)
        {
            names.emplace_back("Thermal Conductivity Dry");
        }
        if(WaterVaporDiffusionResistanceFactor)
        {
            names.emplace_back("Water Vapor Diffusion Resistance Factor");
        }
        if(MoistureStorageFunction)
        {
            names.emplace_back("Moisture Storage Function");
        }
        if(LiquidTransportationSuction)
        {
            names.emplace_back("Liquid Transportation Suction Curve");
        }
        if(LiquidTransportationRedistribution)
        {
            names.emplace_back("Liquid Transportation Redistribution Curve");
        }
        if(ThermalConductivityMoistureAndTemperatureDependent)
        {
            names.emplace_back("Thermal Conductivity Moisture and Temperature Dependent");
        }
        return names;
    }

    std::vector<std::string> MaterialMissingProperties::missingPropertiesMessage() const
    {
        std::vector<std::string> message{};
        const auto names{missingPropertyNames()};
        if(!names.empty())
        {
            message.push_back("Material " + materialName + " is missing following properties:");
            for(const auto & name : names)
            {
                message.push_back("- " + name);
            }
        }
        return message;
    }

    bool MaterialMissingProperties::isMissingAnyProperty() const
    {
        return Density || Emissivity || Porosity || SpecificHeatCapacityDry
               || ThermalConductivityDry || WaterVaporDiffusionResistanceFactor
               || MoistureStorageFunction || LiquidTransportationSuction
               || LiquidTransportationRedistribution
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