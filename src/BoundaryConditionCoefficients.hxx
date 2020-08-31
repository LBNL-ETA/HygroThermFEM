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
    /// VariableBCTARPHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Structure to keep boundary condition coefficients for comprehensive natural
    //! convection model
    //!
    //! It is simply the structure that holds all necessary variables that are needed to create
    //! boundary condition for comprehensive natural convection model that is often named as TARP
    struct VariableBCTARPHCCoefficients
    {
        //! \brief Construction of variables for TARP convection model.
        //! This structure holds only values that vary through every timestep.
        //!
        //! \param airTemperature Air/Ambient temperature of exterior/interior environment
        //! \param airHumidity Air humidity of the environment
        VariableBCTARPHCCoefficients(double airTemperature, double airHumidity);

        double AirTemperature{0.0}; // Celsius
        double AirHumidity{0.0}; // Humidity range is from 0 to 1
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// VariableBCASHRAEInsideHCCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////
    struct VariableBCASHRAEInsideHCCoefficients
    {

        double AirTemperature{0.0}; // Celsius
        double AirHumidity{0.0}; // Humidity range is from 0 to 1
        double SurfaceHeight{0.0}; // meters
        double SurfaceTilt{90.0}; // degrees
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
