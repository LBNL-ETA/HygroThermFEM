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
    static constexpr double STEFANBOLTZMANN = 5.6697E-8;
    static constexpr double ABSOLUTEZERO = -273.15;
    static constexpr double GravityConstant = 9.807;

    // Need to handle phase change through steep function. Melting is happening
    // between IcePoint and FreezingPoint
    static constexpr double FreezingPoint = 0.0;
    static constexpr double IcePoint = -0.001;

    static constexpr double EnthalpyOfFusion = 333550;   // J/kg

    //// TODO: Keep these constant for now. Gases can calculate these properties.
    static constexpr double Density_Air = 1.2922;
    static constexpr double Cp_Air = 1000;
    static constexpr double K_Air = 0.025;

    //// TODO: Check if water and ice properties need to be recalculated
    static constexpr double Density_Water = 1000;
    static constexpr double Cp_Water = 4184;
    static constexpr double K_Water = 0.591;

    static constexpr double Density_Ice = 916.7;
    static constexpr double Cp_Ice = 2108;
    static constexpr double K_Ice = 2.22;

    static constexpr double Cp_Vapor = 1850;

}   // namespace Constants
