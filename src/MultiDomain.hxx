#pragma once

#include "SingleDomain.hxx"

namespace HygroThermFEM
{
    struct MultiDomain
    {
        MultiDomain() = default;
        MultiDomain(bool performThermal, bool performMoisture);

        SingleDomain thermalDomain{DomainType::Thermal};
        SingleDomain moistureDomain{DomainType::Moisture};
        bool simulateThermal{true};
        bool simulateMoisture{true};
    };

}   // namespace HygroThermFEM
