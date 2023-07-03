#pragma once

#include "ThermalDomain.hxx"
#include "MoistureDomain.hxx"

namespace HygroThermFEM
{
    struct MultiDomain
    {
        MultiDomain() = default;
        MultiDomain(bool performThermal, bool performMoisture);

        ThermalDomain thermalDomain;
        MoistureDomain moistureDomain;
        bool simulateThermal{true};
        bool simulateMoisture{true};
    };

}   // namespace HygroThermFEM
