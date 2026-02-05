#include "State.hxx"

namespace HygroThermFEM
{
    double State::getValue(const BaseVariable property) const
    {
        switch(property)
        {
            case BaseVariable::temperature:
                return temperature;
            case BaseVariable::humidity:
                return humidity;
            case BaseVariable::pressure:
                return pressure;
            case BaseVariable::liquidPercent:
                return liquidPercent;
        }
        return 0.0;
    }

    void State::setValue(const BaseVariable property, const double value)
    {
        switch(property)
        {
            case BaseVariable::temperature:
                temperature = value;
                break;
            case BaseVariable::humidity:
                humidity = value;
                break;
            case BaseVariable::pressure:
                pressure = value;
                break;
            case BaseVariable::liquidPercent:
                liquidPercent = value;
                break;
        }
    }

}   // namespace HygroThermFEM
