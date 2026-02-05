#include "State.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////////////////////////
    ///  State
    ////////////////////////////////////////////////////////////////////////////

    State::State() : State(StateParams{})
    {}

    State::State(StateParams params)
    {
        m_Property[BaseVariable::temperature] = params.temperature;
        m_Property[BaseVariable::humidity] = params.humidity;
        m_Property[BaseVariable::pressure] = params.pressure;
        m_Property[BaseVariable::liquidPercent] = params.liquidPercent;
    }

    double State::getValue(const BaseVariable t_Property) const
    {
        return m_Property.at(t_Property);
    }

    void State::setValue(const BaseVariable t_Property, const double t_Value)
    {
        m_Property[t_Property] = t_Value;
    }

}   // namespace HygroThermFEM
