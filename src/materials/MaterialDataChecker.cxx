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

        // Latent heat of fusion is carried by the water a material holds. The water comes from
        // the sorption curve, and porosity splits the vapour off it; a material lacking either
        // is left out of the node water content and would silently freeze nothing. A material
        // that holds no water states so with a zero sorption curve and zero porosity: missing
        // data is not zero.
        const bool fusionNeedsWater{multiDomain.isThermalSimulationON() && !physics.excludeLatentHeatOfFusion};

        // Porosity is indirectly required when calculating liquid content
        missing.Porosity =
          !material.hasPorosity() && (multiDomain.isMoistureSimulationON() || fusionNeedsWater);

        // Used only in case of transient thermal simulation
        missing.SpecificHeatCapacityDry = isTransientSimulation && !material.hasHeatCapacity()
                                          && multiDomain.isThermalSimulationON();

        // A thermal simulation needs a conductivity for every material: its moisture and
        // temperature tables when the option asks for them and the material carries them,
        // its dry value otherwise. So the dry value is missing only when neither is there.
        const bool conductivityFromTables{physics.thermalConductivityMoistureAndTemperatureDependent
                                          && material.hasThermalConductivityMoistureAndTemperatureDependent()};
        missing.ThermalConductivityDry = !material.hasThermalConductivityDry()
                                         && multiDomain.isThermalSimulationON() && !conductivityFromTables;

        // Water vapor diffusion is required in moisture always and in thermal only in case heat of
        // evaporation calculation is on. Either form satisfies it: a material may carry the
        // single value or a curve in water content, and a material carrying only the curve
        // is fully specified.
        missing.WaterVaporDiffusionResistanceFactor =
          (!material.hasAnyDiffusionResistanceFactor() && multiDomain.isMoistureSimulationON())
          || (!material.hasAnyDiffusionResistanceFactor() && !physics.excludeHeatOfEvaporation
              && multiDomain.isThermalSimulationON());

        // Sorption curve: any moisture simulation needs it (the transient capacity term and
        // the water-content evaluation both read it, steady included), and a thermal
        // simulation needs it when capillary conduction or latent heat of fusion is on.
        missing.MoistureStorageFunction =
          (!material.hasSorptionCurve() && multiDomain.isThermalSimulationON()
           && !physics.excludeCapillaryConduction)
          || (!material.hasSorptionCurve() && multiDomain.isMoistureSimulationON())
          || (!material.hasSorptionCurve() && fusionNeedsWater);


        // Two options in case of Liquid Transportation Curve
        //  1. Thermal simulation includes capillary conduction.
        //  2. Moisture simulation include capillary transportation.
        missing.LiquidTransportationSuction =
          (!material.hasLiquidTransportationCurve() && multiDomain.isThermalSimulationON()
           && !physics.excludeCapillaryConduction)
          || (!material.hasLiquidTransportationCurve() && multiDomain.isMoistureSimulationON()
              && !physics.excludeWaterLiquidTransportation);

        // missing.LiquidTransportationRedistribution =

        // The tables are never required on their own: a material without them falls back
        // to its dry conductivity, and a material with neither is reported above.
        missing.ThermalConductivityMoistureAndTemperatureDependent = false;

        return missing;
    }
}   // namespace HygroThermFEM