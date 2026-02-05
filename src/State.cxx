#include "State.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///  State
    ////////////////////////////////////////////////////////////////////////////

    State::State(const double t_Temperature,
                 const double t_Humidity,
                 const double t_Pressure,
                 const double liquidPercent)
    {
        m_Property[BaseVariable::temperature] = t_Temperature;
        m_Property[BaseVariable::humidity] = t_Humidity;
        m_Property[BaseVariable::pressure] = t_Pressure;
        m_Property[BaseVariable::liquidPercent] = liquidPercent;
    }

    State::State(StateParams params)
        : State(params.temperature, params.humidity, params.pressure, params.liquidPercent)
    {}

    double State::getValue(const BaseVariable t_Property) const
    {
        return m_Property.at(t_Property);
    }

    void State::setValue(const BaseVariable t_Property, const double t_Value)
    {
        m_Property[t_Property] = t_Value;
    }

}   // namespace HygroThermFEM
