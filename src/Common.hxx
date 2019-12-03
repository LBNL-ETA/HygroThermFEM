#pragma once

#include <cmath>

namespace HygroThermFEM
{
    //! \brief Conversion of temperature from Celsius to Kelvin
    //!
    //! \param temperature Input temperature in degrees Celsius
    inline double toKelvin(const double temperature)
    {
        return temperature + 273.15;
    }

}   // namespace HygroThermFEM

namespace Constants
{
    static const double PI_HTF = atan(1) * 4;
    static const double STEFANBOLTZMANN = 5.6697E-8;
    static const double ABSOLUTEZERO = -273.15;

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

namespace TimestepData
{
    //! Enumeration that will be used to notify exterior world on current timestep simlation level.
    //!
    //! Standard - current timestep level is normal. No division performed and engine is trying to do
    //! simulation with provided timestep division
    //! One - current timestep level is increased. This case means that the engine is trying to
    //! perform timestep simulation by dividing it with some number (current implementation divides
    //! timestep into ten smaller timesteps) Two - Timestep is even more refined. Three - Maximum
    //! level of timestep refinement.

    enum class Level
    {
        Standard = 0,
        One,
        Two,
        Three
    };

    inline Level operator++(Level& mode)
    {
        static unsigned numberOfEnums{static_cast<int>(Level::Three) + 1};
        mode = static_cast<Level>((static_cast<int>(mode) + 1) % numberOfEnums);
        return mode;
    }

    //! Engine will try to divide timestep three times before it reports divergence.
    const unsigned maxDivisions{3};

    //! In case of convergence failure, program will divide timestep into number of smaller timesteps.
    const unsigned NumberOfSubsegments{10};
}   // namespace Timestep
