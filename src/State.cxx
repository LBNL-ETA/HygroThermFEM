#include "State.hxx"

namespace MoisThermFEM
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

    State operator+(const State & lhs, const State & rhs)
    {
        State state(0, 0, 0, 0);
        state.setValue(BaseVariable::temperature,
                       lhs.getValue(BaseVariable::temperature) + rhs.getValue(BaseVariable::temperature));
        state.setValue(BaseVariable::humidity,
                       lhs.getValue(BaseVariable::humidity) + rhs.getValue(BaseVariable::humidity));
        state.setValue(BaseVariable::pressure,
                       lhs.getValue(BaseVariable::pressure) + rhs.getValue(BaseVariable::pressure));
        state.setValue(BaseVariable::liquidPercent,
                       lhs.getValue(BaseVariable::liquidPercent)
                         + rhs.getValue(BaseVariable::liquidPercent));

        return state;
    }

    State operator-(const State & lhs, const State & rhs)
    {
        State state(0, 0, 0, 0);
        state.setValue(BaseVariable::temperature,
                       lhs.getValue(BaseVariable::temperature) - rhs.getValue(BaseVariable::temperature));
        state.setValue(BaseVariable::humidity,
                       lhs.getValue(BaseVariable::humidity) - rhs.getValue(BaseVariable::humidity));
        state.setValue(BaseVariable::pressure,
                       lhs.getValue(BaseVariable::pressure) - rhs.getValue(BaseVariable::pressure));
        state.setValue(BaseVariable::liquidPercent,
                       lhs.getValue(BaseVariable::liquidPercent)
                         - rhs.getValue(BaseVariable::liquidPercent));
        return state;
    }

    double State::getValue(const BaseVariable t_Property) const
    {
        return m_Property.at(t_Property);
    }

    void State::setValue(const BaseVariable t_Property, const double t_Value)
    {
        m_Property[t_Property] = t_Value;
    }

}   // namespace MoisThermFEM
