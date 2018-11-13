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
        m_Property[StateProperty::temperature] = t_Temperature;
        m_Property[StateProperty::humidity] = t_Humidity;
        m_Property[StateProperty::pressure] = t_Pressure;
        m_Property[StateProperty::liquidPercent] = liquidPercent;
    }

    State operator+(const State & lhs, const State & rhs)
    {
        State state(0, 0, 0, 0);
        state.setValue(StateProperty::temperature,
                       lhs.getValue(StateProperty::temperature) + rhs.getValue(StateProperty::temperature));
        state.setValue(StateProperty::humidity,
                       lhs.getValue(StateProperty::humidity) + rhs.getValue(StateProperty::humidity));
        state.setValue(StateProperty::pressure,
                       lhs.getValue(StateProperty::pressure) + rhs.getValue(StateProperty::pressure));
        state.setValue(StateProperty::liquidPercent,
                       lhs.getValue(StateProperty::liquidPercent)
                         + rhs.getValue(StateProperty::liquidPercent));

        return state;
    }

    State operator-(const State & lhs, const State & rhs)
    {
        State state(0, 0, 0, 0);
        state.setValue(StateProperty::temperature,
                       lhs.getValue(StateProperty::temperature) - rhs.getValue(StateProperty::temperature));
        state.setValue(StateProperty::humidity,
                       lhs.getValue(StateProperty::humidity) - rhs.getValue(StateProperty::humidity));
        state.setValue(StateProperty::pressure,
                       lhs.getValue(StateProperty::pressure) - rhs.getValue(StateProperty::pressure));
        state.setValue(StateProperty::liquidPercent,
                       lhs.getValue(StateProperty::liquidPercent)
                         - rhs.getValue(StateProperty::liquidPercent));
        return state;
    }

    double State::getValue(const StateProperty t_Property) const
    {
        return m_Property.at(t_Property);
    }

    void State::setValue(const StateProperty t_Property, const double t_Value)
    {
        m_Property[t_Property] = t_Value;
    }

    double State::getLiquidPercent() const
    {
        return m_Property.at(StateProperty::liquidPercent);
    }

}   // namespace MoisThermFEM
