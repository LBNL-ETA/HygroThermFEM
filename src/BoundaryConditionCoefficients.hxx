#pragma once

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // FixedBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

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

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // VariableBCHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Structure to keep boundary condition coefficients for variable heat transfer
    //! coefficients case
    struct VariableBCHCCoefficients
    {
        VariableBCHCCoefficients(double airTemperature, double airHumidity);

        double AirTemperature{0};
        double AirHumidity{0};
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ConstantBCTemperatures
    ///////////////////////////////////////////////////////////////////////////////////////////////

    struct ConstantBCTemperatures
    {
        ConstantBCTemperatures(double temperature1, double temperature2);

        double Temperature1{0};
        double Temperature2{0};
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // LinearizedRadiationBCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    struct LinearizedRadiationBCCoefficients
    {
        LinearizedRadiationBCCoefficients(double radiationCoefficient, double radiationTemperature);

        double RadiationCoefficient{0};
        double RadiationTemperature{0};
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // BlackBodyRadiationBCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    struct BlackBodyRadiationBCCoefficients
    {
        BlackBodyRadiationBCCoefficients(double emissivity, double temperature);

        double Emissivity{0.0};
        double Temperature{0.0};
    };

}   // namespace HygroThermFEM
