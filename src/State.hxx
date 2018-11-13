#pragma once

#include <map>

namespace MoisThermFEM
{
    enum class StateProperty
    {
        temperature,
        humidity,
        pressure,
        liquidPercent
    };

    // Simple class to hold state variables
    class State
    {
    public:
        explicit State(const double t_Temperature = 0,
                       const double t_Humidity = 0,
                       const double t_Pressure = 101325,
                       const double liquidPercent = 1.0);

        friend State operator+(const State & lhs, const State & rhs);
        friend State operator-(const State & lhs, const State & rhs);

        double getValue(const StateProperty t_Property) const;
        void setValue(const StateProperty t_Property, const double t_Value);
        double getLiquidPercent() const;

    private:
        std::map<StateProperty, double> m_Property;
    };

}   // namespace MoisThermFEM