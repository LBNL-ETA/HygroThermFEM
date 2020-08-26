#pragma once

#include <cmath>

namespace HygroThermFEM
{
    //! \brief Conversion of temperature from Celsius to Kelvin
    //!
    //! \param temperature Input temperature in degrees Celsius
    inline double celsiusToKelvin(const double temperature)
    {
        return temperature + 273.15;
    }

    enum class SimulationType{SteadyState, Transient};

}   // namespace HygroThermFEM

namespace Constants
{
    static const double PI_HTF = atan(1) * 4;
    static const double STEFANBOLTZMANN = 5.6697E-8;
    static const double ABSOLUTEZERO = -273.15;
    static const double GravityConstant = 9.807;

    // Need to handle phase change through steep function. Melting is happening
    // between IcePoint and FreezingPoint
    static const double FreezingPoint = 0.0;
    static const double IcePoint = -0.001;

    static const double EnthalpyOfFusion = 333550;   // J/kg

    //// TODO: Keep these constant for now. Gases can calculate these properties.
    static const double Density_Air = 1.2922;
    static const double Cp_Air = 1000;
    static const double K_Air = 0.025;

    //// TODO: Check if water and ice properties need to be recalculated
    static const double Density_Water = 1000;
    static const double Cp_Water = 4184;
    static const double K_Water = 0.591;

    static const double Density_Ice = 916.7;
    static const double Cp_Ice = 2108;
    static const double K_Ice = 2.22;

    static const double Cp_Vapor = 1850;

}   // namespace Constants
