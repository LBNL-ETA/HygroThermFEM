#include "MaterialDataChecker.hxx"
#include "MultiDomain.hxx"
#include "MaterialPool.hxx"
#include "Material.hxx"
#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    MaterialDataChecker::MaterialDataChecker(const MultiDomain & multiDomain) :
        multiDomain(multiDomain)
    {}

    std::vector<MaterialMissingProperties>
      MaterialDataChecker::checkMaterialProperties(const bool isTransientSimulation)
    {
        std::vector<MaterialMissingProperties> missingProperties;
        for(const auto & materialName : MaterialPool::Instance().getSolidMaterials())
        {
            const auto & material{MaterialPool::Instance().material(materialName)};
            const auto materialCheck = checkMaterial(material, isTransientSimulation);
            if(materialCheck.isMissingAnyProperty())
            {
                missingProperties.push_back(materialCheck);
            }
        }

        return missingProperties;
    }

    MaterialMissingProperties MaterialDataChecker::checkMaterial(const IMaterial & material,
                                                                 const bool isTransientSimulation)
    {
        MaterialMissingProperties missing;
        missing.Density =
          isTransientSimulation && !material.hasDensity() && multiDomain.isMoistureSimulationON();
        missing.Emissivity = !material.hasEmissivity();
        missing.Porosity =
          isTransientSimulation && !material.hasPorosity() && multiDomain.isMoistureSimulationON();
        missing.SpecificHeatCapacityDry = isTransientSimulation && !material.hasHeatCapacity();
        missing.ThermalConductivityDry =
          !material.hasThermalConductivityDry() && multiDomain.isThermalSimulationON()
          && !SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent();
        missing.WaterVaporDiffusionResistanceFactor =
          !material.hasDiffusionResistanceFactor() && multiDomain.isMoistureSimulationON();
        missing.MoistureStorageFunction = isTransientSimulation && !material.hasSorptionCurve()
                                          && multiDomain.isMoistureSimulationON();
        missing.LiquidTransportationSuction =
          !material.hasLiquidTransportationCurve() && multiDomain.isMoistureSimulationON();
        // missing.LiquidTransportationRedistribution =
        missing.ThermalConductivityMoistureAndTemperatureDependent =
          !material.hasThermalConductivityMoistureAndTemperatureDependent()
          && multiDomain.isThermalSimulationON();
        return missing;
    }

}   // namespace HygroThermFEM