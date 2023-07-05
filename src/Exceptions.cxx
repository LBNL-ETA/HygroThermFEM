#include "Exceptions.hxx"

namespace HygroThermFEM
{
    InvalidSimulationTypeException::InvalidSimulationTypeException() :
        std::runtime_error("Invalid simulation type selected.")
    {}

    SolutionFailedToConvergeException::SolutionFailedToConvergeException() :
        std::runtime_error("Solution failed to converge.")
    {}
}   // namespace HygroThermFEM