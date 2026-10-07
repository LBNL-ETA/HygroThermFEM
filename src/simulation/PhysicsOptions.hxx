#pragma once

namespace HygroThermFEM
{
    //! \brief Immutable snapshot of the physics-model flags used by a single domain solve.
    //!
    //! These values were historically read from the process-global SimulationProperties
    //! singleton at element construction and assembly time. Passing them as an explicit value
    //! allows per-domain physics and isolates tests from shared global state. When a domain is
    //! not given an explicit instance it falls back to fromGlobals(), so the default behaviour
    //! is unchanged.
    //!
    //! Element constructors read these flags to decide which equation terms to register, so a
    //! domain must receive its options BEFORE its geometry is built (setPhysicsOptions before
    //! the first createElement call).
    struct PhysicsOptions
    {
        //! Excludes water liquid transportation from the mass transfer equation.
        bool excludeWaterLiquidTransportation{false};
        //! Excludes heat of evaporation from the heat transfer equation.
        bool excludeHeatOfEvaporation{false};
        //! Excludes capillary conduction from the heat transfer equation.
        bool excludeCapillaryConduction{false};
        //! Excludes water vapor diffusion conduction from the heat transfer equation.
        bool excludeVaporDiffusionConduction{false};
        //! Includes moisture and temperature dependency of thermal conductivity.
        bool thermalConductivityMoistureAndTemperatureDependent{false};
        //! Excludes the latent-heat-of-fusion (freezing) capacity from the heat equation.
        //! Excluded by default, matching the engine's historical behaviour.
        bool excludeLatentHeatOfFusion{true};

        //! \brief Builds the options from the process-global singleton (the historical source).
        [[nodiscard]] static PhysicsOptions fromGlobals();
    };
}   // namespace HygroThermFEM
