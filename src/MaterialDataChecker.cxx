#include "MaterialDataChecker.hxx"
#include "MultiDomain.hxx"
#include "MaterialPool.hxx"
#include "Material.hxx"
#include "SimulationProperties.hxx"
#include "Exceptions.hxx"

namespace HygroThermFEM
{

    namespace
    {
        //! \brief Check single material properties against current engine settings
        //!
        //! \param material Material that needs to be checked.
        //! \param isTransientSimulation True if simulation is transient, False if steady-state.
        //! \return Missing properties for given material.
        MaterialMissingProperties checkMaterial(const MultiDomain & domain,
                                                const IMaterial & material,
                                                bool isTransientSimulation)
        {
            MaterialMissingProperties missing;

            missing.materialName = material.name();

            // Used in thermal equation only in case of transient simulation
            missing.Density = isTransientSimulation && !material.hasDensity() && domain.simulateThermal;

            // Emissivity will be required always
            missing.Emissivity = !material.hasEmissivity() && domain.simulateThermal;

            // Porosity is indirectly required when calculating liquid content
            missing.Porosity = !material.hasPorosity() && domain.simulateMoisture;

            // Used only in case of transient thermal simulation
            missing.SpecificHeatCapacityDry =
              isTransientSimulation && !material.hasHeatCapacity() && domain.simulateThermal;

            // Thermal conductivity dry (single value) is required only in thermal simulation. Also, in
            // case when thermal conductivity moisture and temperature dependent is used, then this is
            // not required property.
            missing.ThermalConductivityDry =
              !material.hasThermalConductivityDry() && domain.simulateThermal
              && !SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent();

            // Water vapor diffusion is required in moisture always and in thermal only in case heat of
            // evaporation calculation is on
            missing.WaterVaporDiffusionResistanceFactor =
              (!material.hasDiffusionResistanceFactor() && domain.simulateMoisture)
              || (!material.hasDiffusionResistanceFactor()
                  && !SimulationProperties::Instance().excludeHeatOfEvaporation()
                  && domain.simulateThermal);

            // Three options in case of Sorption Curve:
            //  1. Thermal simulation includes capillary conduction.
            //  2. Moisture simulation include capillary transportation.
            //  3. Moisture simulation is transient.
            missing.MoistureStorageFunction =
              (!material.hasSorptionCurve() && domain.simulateThermal
               && !SimulationProperties::Instance().excludeCapillaryConduction())
              || (!material.hasSorptionCurve() && domain.simulateMoisture
                  && !SimulationProperties::Instance().excludeWaterLiquidTransportation())
              || (!material.hasSorptionCurve() && isTransientSimulation && domain.simulateMoisture);


            // Two options in case of Liquid Transportation Curve
            //  1. Thermal simulation includes capillary conduction.
            //  2. Moisture simulation include capillary transportation.
            missing.LiquidTransportationSuction =
              (!material.hasLiquidTransportationCurve() && domain.simulateThermal
               && !SimulationProperties::Instance().excludeCapillaryConduction())
              || (!material.hasLiquidTransportationCurve() && domain.simulateMoisture
                  && !SimulationProperties::Instance().excludeWaterLiquidTransportation());

            // missing.LiquidTransportationRedistribution =

            missing.ThermalConductivityMoistureAndTemperatureDependent =
              !material.hasThermalConductivityMoistureAndTemperatureDependent()
              && domain.simulateThermal
              && SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent();

            return missing;
        }

        //! \brief Check all material properties against current engine settings.
        //!
        //! \param isTransientSimulation Different checks are needed for different simulation types.
        //! \return All missing properties for every material.
        MaterialsErrorCheckVector checkMaterialProperties(const MultiDomain & domain,
                                                          const bool isTransientSimulation)
        {
            MaterialsErrorCheckVector missingProperties;
            for(const auto & materialName : MaterialPool::Instance().getSolidMaterials())
            {
                const auto & material{MaterialPool::Instance().material(materialName)};
                const auto materialCheck = checkMaterial(domain, material, isTransientSimulation);
                if(materialCheck.isMissingAnyProperty())
                {
                    missingProperties.push_back(materialCheck);
                }
            }

            return missingProperties;
        }

        MaterialsErrorCheckVector checkMaterialsForTransientSimulation(const MultiDomain & domain)
        {
            return checkMaterialProperties(domain, true);
        }

        MaterialsErrorCheckVector checkMaterialsForSteadyStateSimulation(const MultiDomain & domain)
        {
            return checkMaterialProperties(domain, false);
        }

        const std::map<SimulationType,
                       std::function<MaterialsErrorCheckVector(const MultiDomain &)>>
          simulationCheckers{
            {SimulationType::SteadyState, checkMaterialsForSteadyStateSimulation},
            {SimulationType::Transient, checkMaterialsForTransientSimulation},
          };
    }   // namespace

    MaterialsErrorCheckVector checkForMaterialsValidity(const MultiDomain & domain,
                                                        const SimulationType simulationType)
    {
        auto it = simulationCheckers.find(simulationType);
        if(it == simulationCheckers.end())
        {
            throw InvalidSimulationTypeException();
        }
        return it->second(domain);
    }

}   // namespace HygroThermFEM