#pragma once

namespace HygroThermFEM
{
    //! \brief Structure to keep boundary condition coefficients for fixed heat transfer
    //! coefficients case
    struct FixedBCHCCoefficients
    {
        FixedBCHCCoefficients(double airTemperature, double convectionCoefficient);

        FixedBCHCCoefficients(double airTemperature,
                              double convectionCoefficient,
                              double airHumidity);

        double AirTemperature{0};
        double ConvectionCoefficient{0};
        double AirHumidity{0};
    };

    //! \brief Structure to keep boundary condition coefficients for variable heat transfer
    //! coefficients case
    struct VariableBCHCCoefficients
    {
        VariableBCHCCoefficients(double airTemperature, double airHumidity);

        double AirTemperature{0};
        double AirHumidity{0};
    };
}   // namespace HygroThermFEM
