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

    // Phase change happens over the ramp [IcePoint, FreezingPoint]: the liquid fraction
    // goes linearly from 0 to 1 across it. Physically, pore water freezes over a range
    // (freezing-point depression), not at a single temperature. Numerically the width
    // matters a great deal: with a 1 mK ramp the fusion capacity is a 4-orders-of-
    // magnitude cliff that PINS crossing nodes at its edge and 2-cycles the NR iteration
    // at every freeze-thaw crossing (each falsely-accepted cycle midpoint leaks energy);
    // at 0.1 K the mushy band is a resolvable, continuous region the iteration can
    // converge into.
    static constexpr double FreezingPoint = 0.0;
    static constexpr double IcePoint = -0.1;

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
