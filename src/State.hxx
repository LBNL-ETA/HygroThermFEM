#pragma once

namespace HygroThermFEM
{
    //! \brief Enumerator that holds basic state variables defined in finite element model
    //!
    //! Problem is solved for three base state variables that are solved independently through
    //! iterations.
    enum class BaseVariable
    {
        temperature,    //!< Temperature
        humidity,       //!< Humidity
        pressure,       //!< Pressure
        liquidPercent   //!< Percent of water content that is in liquid state
    };

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

        //! Returns state value for given BaseVariable (for runtime access)
        [[nodiscard]] double getValue(BaseVariable property) const;

        //! Sets state value for given BaseVariable (for runtime access)
        void setValue(BaseVariable property, double value);
    };

}   // namespace HygroThermFEM
