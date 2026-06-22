#include "SolverSettings.hxx"

#include "SimulationProperties.hxx"
#include "TimestepData.hxx"

namespace HygroThermFEM
{
    SolverSettings SolverSettings::fromGlobals()
    {
        const auto & properties = SimulationProperties::Instance();
        const auto & timesteps = Timesteps::Settings::Instance();

        return SolverSettings{.relaxationParameter = properties.relaxationParamter(),
                              .errorTolerance = properties.errorTolerance(),
                              .maxNumberOfIterations = properties.maxNumberOfIterations(),
                              .maxDivisions = timesteps.getMaxDivisions(),
                              .numberOfSubtimesteps = timesteps.getNumberOfSubtimesteps()};
    }
}   // namespace HygroThermFEM
