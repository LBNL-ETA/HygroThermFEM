#include "PhysicsOptions.hxx"

#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    PhysicsOptions PhysicsOptions::fromGlobals()
    {
        const auto & properties = SimulationProperties::Instance();

        return PhysicsOptions{
          .excludeWaterLiquidTransportation = properties.excludeWaterLiquidTransportation(),
          .excludeHeatOfEvaporation = properties.excludeHeatOfEvaporation(),
          .excludeCapillaryConduction = properties.excludeCapillaryConduction(),
          .excludeVaporDiffusionConduction = properties.excludeVaporDiffusionConduction(),
          .thermalConductivityMoistureAndTemperatureDependent =
            properties.thermalConductivityTemperatureAndMoistureDependent(),
          .excludeLatentHeatOfFusion = properties.excludeLatentHeatOfFusion()};
    }
}   // namespace HygroThermFEM
