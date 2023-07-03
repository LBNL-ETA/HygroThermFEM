#pragma once

#include <stdexcept>

namespace HygroThermFEM
{
    class InvalidSimulationTypeException : public std::runtime_error
    {
    public:
        InvalidSimulationTypeException();
    };
}