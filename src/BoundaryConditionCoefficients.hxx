#pragma once

namespace HygroThermFEM
{
    //! \brief Structure to keep boundary condition coefficients for fixed heat transfer
    //! coefficients case
    struct FixedBCHCCoefficients
    {
        FixedBCHCCoefficients(const double airTemperature, const double convectionCoefficient) :
            AirTemperature(airTemperature), ConvectionCoefficient(convectionCoefficient)
        {}

        FixedBCHCCoefficients(const double airTemperature,
                              const double convectionCoefficient,
                              const double airHumidity) :
            AirTemperature(airTemperature),
            ConvectionCoefficient(convectionCoefficient),
            AirHumidity(airHumidity)
        {}

        double AirTemperature{0};
        double ConvectionCoefficient{0};
        double AirHumidity{0};
    };
}   // namespace HygroThermFEM
