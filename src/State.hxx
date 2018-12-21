#pragma once

#include <map>

namespace MoisThermFEM
{
    //! \brief Enumerator that hold basic state variables defined in finite element model
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
    //! Basic state variables are temperature, humidity and pressure. In addition, class hold
    //! important information about percentage of water content that is in liquid state.
    class State
    {
    public:
        //! Construction of State.
        explicit State(
          double t_Temperature = 0,     //!< Temperature (in Celsius).
          double t_Humidity = 0,        //!< Humidity (range from zero to one).
          double t_Pressure = 101325,   //!< Pressure.
          double liquidPercent = 1.0    //!< Percent of water content in liquid state.
        );

        //! Operator overloading for addition operation
        friend State operator+(const State & lhs, const State & rhs);

        //! Operator overloading for subtraction operation
        friend State operator-(const State & lhs, const State & rhs);

        //! Returns state value for given BaseVariable
        double getValue(BaseVariable t_Property) const;

        //! Sets state value for given BaseVariable
        void setValue(
          BaseVariable t_Property,   //!< BaseVariable for which value will be set to.
          double t_Value              //!< New value of given BaseVariable
        );

    private:
        std::map<BaseVariable, double> m_Property;
    };

}   // namespace MoisThermFEM