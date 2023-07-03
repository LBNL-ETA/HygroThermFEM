#include "Exceptions.hxx"

namespace HygroThermFEM
{
    InvalidSimulationTypeException::InvalidSimulationTypeException() :
        std::runtime_error("Invalid simulation type selected.")
    {}
}   // namespace HygroThermFEM