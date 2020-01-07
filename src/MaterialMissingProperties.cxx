#include "MaterialMissingProperties.hxx"

namespace HygroThermFEM
{
    //////////////////////////////////////////////////////////////////////////
    ///  MaterialMissingProperties
    //////////////////////////////////////////////////////////////////////////
    std::string MaterialMissingProperties::missingPropertiesMessage() const
    {
        std::string message{};
        if(isMissingAnyProperty())
        {
            message += "Material " + materialName + " is missing following properties:\n";
            if(Density)
            {
                message += "- Density \n";
            }
            if(Emissivity)
            {
                message += "- Emissivity \n";
            }
            if(Porosity)
            {
                message += "- Porosity \n";
            }
            if(SpecificHeatCapacityDry)
            {
                message += "- Specific Heat Capacity Dry \n";
            }
            if(ThermalConductivityDry)
            {
                message += "- Thermal Conductivity Dry \n";
            }
            if(WaterVaporDiffusionResistanceFactor)
            {
                message += "- Water Vapor Diffusion Resistance Factor \n";
            }
            if(MoistureStorageFunction)
            {
                message += "- Moisture Storage Function \n";
            }
            if(LiquidTransportationSuction)
            {
                message += "- Liquid Transportation Suction Curve \n";
            }
            if(LiquidTransportationRedistribution)
            {
                message += "- Liquid Transportation Redistribution Curve \n";
            }
            if(ThermalConductivityMoistureAndTemperatureDependent)
            {
                message += "- Thermal Conductivity Moisture and Temperature Dependent \n";
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