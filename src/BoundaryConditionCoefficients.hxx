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
    /// TARPCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Structure to keep boundary condition coefficients for comprehensive natural
    //! convection model
    //!
    //! It is simply the structure that holds all necessary variables that are needed to create
    //! boundary condition for comprehensive natural convection model that is often named as TARP
    struct TARPCoefficients
    {
        //! \brief Construction of variables for TARP convection model.
        //! This structure holds only values that vary through every timestep.
        //!
        //! \param airTemperature Air/Ambient temperature of exterior/interior environment
        //! \param airHumidity Air humidity of the environment
        TARPCoefficients(double airTemperature, double airHumidity);

        double AirTemperature{0.0};   // Celsius
        double AirHumidity{0.0};      // Humidity range is from 0 to 1
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// ASHRAEInsideCoefficients
    ///////////////////////////////////////////////////////////////////////////////////////////////

    //! \brief Structure to keep boundary condition coefficients that are used in ASHRAE indoor
    //! model.
    //!
    //! This structure will keep only coefficients that are variable between timesteps. Idea is that
    //! this structure will be used to keep data in some kind of array.
    struct ASHRAEInsideCoefficients
    {
        //! \brief Simple constructor for the structure
        ASHRAEInsideCoefficients(double air_temperature, double air_pressure, double air_humidity);

        double AirTemperature{0.0};     // Celsius
        double AirPressure{101325.0};   // Pascals
        double AirHumidity{0.0};        // Humidity range is from 0 to 1
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
