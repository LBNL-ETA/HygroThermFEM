#pragma once

namespace HygroThermFEM
{
    //! \brief Holds state variables that finite element model is solved for.
    //!
    //! Basic state variables are temperature, humidity and pressure. In addition, struct holds
    //! important information about percentage of water content that is in liquid state.
    struct State
    {
        double temperature = 0.0;
        double humidity = 0.0;
        double pressure = 101325.0;
        double liquidPercent = 1.0;
    };

}   // namespace HygroThermFEM
