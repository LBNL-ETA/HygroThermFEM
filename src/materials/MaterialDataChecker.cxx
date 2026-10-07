#include "MaterialDataChecker.hxx"
#include "MultiDomain.hxx"
#include "Materials.hxx"
#include "Material.hxx"

namespace HygroThermFEM
{
    //////////////////////////////////////////////////////////////////////////
    ///  MaterialDataChecker
    //////////////////////////////////////////////////////////////////////////

    MaterialDataChecker::MaterialDataChecker(const MultiDomain & multiDomain) :
        multiDomain(multiDomain)
    {}

    MaterialsErrorCheckVector
      MaterialDataChecker::checkMaterialProperties(const bool isTransientSimulation)
    {
        MaterialsErrorCheckVector missingProperties;
        const auto & materialPool = multiDomain.materials();
        for(const auto & materialName : materialPool.getSolidMaterials())
        {
            const auto & material{materialPool.material(materialName)};
            const auto materialCheck = checkMaterial(material, isTransientSimulation);
            if(materialCheck.isMissingAnyProperty())
            {
                missingProperties.push_back(materialCheck);
            }
        }

        return missingProperties;
    }

    MaterialMissingProperties MaterialDataChecker::checkMaterial(const IMaterial & material,
                                                                 const bool isTransientSimulation) const
    {
        MaterialMissingProperties missing;

        const auto physics{multiDomain.physicsOptions()};

        missing.materialName = material.name();

        // Used in thermal equation only in case of transient simulation
        missing.Density =
          isTransientSimulation && !material.hasDensity() && multiDomain.isThermalSimulationON();

        // Emissivity will be required always
        missing.Emissivity = !material.hasEmissivity() && multiDomain.isThermalSimulationON();

        // Porosity is indirectly required when calculating liquid content
        missing.Porosity = !material.hasPorosity() && multiDomain.isMoistureSimulationON();

        // Used only in case of transient thermal simulation
        missing.SpecificHeatCapacityDry = isTransientSimulation && !material.hasHeatCapacity()
                                          && multiDomain.isThermalSimulationON();

        // Thermal conductivity dry (single value) is required only in thermal simulation. Also, in
        // case when thermal conductivity moisture and temperature dependent is used, then this is
        // not required property.
        missing.ThermalConductivityDry =
          !material.hasThermalConductivityDry() && multiDomain.isThermalSimulationON()
          && !physics.thermalConductivityMoistureAndTemperatureDependent;

        // Water vapor diffusion is required in moisture always and in thermal only in case heat of
        // evaporation calculation is on
        missing.WaterVaporDiffusionResistanceFactor =
          (!material.hasDiffusionResistanceFactor() && multiDomain.isMoistureSimulationON())
          || (!material.hasDiffusionResistanceFactor() && !physics.excludeHeatOfEvaporation
              && multiDomain.isThermalSimulationON());

        // Sorption curve: any moisture simulation needs it (the transient capacity term and
        // the water-content evaluation both read it, steady included), and a thermal
        // simulation needs it when capillary conduction is on.
        missing.MoistureStorageFunction =
          (!material.hasSorptionCurve() && multiDomain.isThermalSimulationON()
           && !physics.excludeCapillaryConduction)
          || (!material.hasSorptionCurve() && multiDomain.isMoistureSimulationON());


        // Two options in case of Liquid Transportation Curve
        //  1. Thermal simulation includes capillary conduction.
        //  2. Moisture simulation include capillary transportation.
        missing.LiquidTransportationSuction =
          (!material.hasLiquidTransportationCurve() && multiDomain.isThermalSimulationON()
           && !physics.excludeCapillaryConduction)
          || (!material.hasLiquidTransportationCurve() && multiDomain.isMoistureSimulationON()
              && !physics.excludeWaterLiquidTransportation);

        // missing.LiquidTransportationRedistribution =

        missing.ThermalConductivityMoistureAndTemperatureDependent =
          !material.hasThermalConductivityMoistureAndTemperatureDependent()
          && multiDomain.isThermalSimulationON()
          && physics.thermalConductivityMoistureAndTemperatureDependent;

        return missing;
    }
}   // namespace HygroThermFEM