#pragma once

#include <map>

namespace MoisThermFEM
{
    enum class Property
    {
        temperature,
        humidity,
        pressure,
        liquidPercent
    };

    enum class Iteration
    {
        Current,
        Previous
    };

    // Simple class to hold state variables
    class State
    {
    public:
        explicit State(const double t_Temperature = 0,
                       const double t_Humidity = 0,
                       const double t_Pressure = 101325,
                       const double liquidPercent = 1.0);

        State(const State & other) = default;

        State & operator=(const State & other) = default;

        friend State operator+(const State & lhs, const State & rhs);
        friend State operator-(const State & lhs, const State & rhs);

        double getValue(const Property t_Property,
                        const Iteration t_Iteration = Iteration::Current) const;
        void setValue(const Property t_Property, const double t_Value);
        double getDeltaValue(const Property t_Property) const;
        double getLiquidPercent() const;

        const std::map<Property, double> & getCurrentValues() const;

    private:
        std::map<Iteration, std::map<Property, double>> m_Property;

        /// Keeps fraction of liquid water that is in liquid state. Opposite would
        /// be frozen state.
        double m_LiquidPercent;
    };

}   // namespace MoisThermFEM