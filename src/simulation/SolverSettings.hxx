#pragma once

#include <cstddef>

namespace HygroThermFEM
{
    //! \brief Immutable snapshot of the solver configuration used by a single domain solve.
    //!
    //! These values were historically read from the process-global SimulationProperties and
    //! Timesteps::Settings singletons inside the solve loop. Passing them as an explicit value
    //! allows per-domain configuration and isolates tests from shared global state. When a domain
    //! is not given an explicit instance it falls back to fromGlobals(), so the default behaviour
    //! is unchanged.
    struct SolverSettings
    {
        double relaxationParameter;
        double errorTolerance;
        std::size_t maxNumberOfIterations;
        unsigned maxDivisions;
        unsigned numberOfSubtimesteps;

        //! \brief Builds the settings from the process-global singletons (the historical source).
        [[nodiscard]] static SolverSettings fromGlobals();
    };
}   // namespace HygroThermFEM
