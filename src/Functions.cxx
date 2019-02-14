#include <algorithm>

#include "Functions.hxx"
#include "Node2D.hxx"
#include "Common.hxx"

namespace HygroThermFEM
{
    double vaporPressureAtTemperature(const double t_temperature)
    {
        const auto temperature = t_temperature + 273.15;
        return std::exp(77.345 + 0.0057 * temperature - 7235.0 / temperature)
               / std::pow(temperature, 8.2);
    }

    double saturationConcentrationAtTemperature(const double t_temperature)
    {
        // RT/M for water vapor
        const auto gasConstantForWaterVapor = 461.4;
        return vaporPressureAtTemperature(t_temperature)
               / ((t_temperature + 273.15) * gasConstantForWaterVapor);
    }

    //////////////////////////////////////////////////////////////////
    ///  IValue
    //////////////////////////////////////////////////////////////////

    std::vector<double> IValue::values(const INodes & nodes) const
    {
        std::vector<double> result(nodes.size(), 0);
        for(size_t i = 0u; i < nodes.size(); ++i)
        {
            result[i] = value(nodes[i]);
        }
        return result;
    }

    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    IFunction::IFunction(const Variable t_Property) : m_Property(t_Property)
    {}

    double IFunction::value(const Node2D & node) const
    {
        return evaluateFunction(node.property(m_Property),
                                node.property(m_Property, Timestep::Previous));
    }

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    Constant::Constant(const double value) : IFunction(Variable::temperature), m_Value(value)
    {}

    double Constant::evaluateFunction(const double, const double) const
    {
        return m_Value;
    }

    //////////////////////////////////////////////////////////////////
    ///  State value
    //////////////////////////////////////////////////////////////////

    StateValue::StateValue(Variable property) : IFunction(property)
    {}

    double StateValue::evaluateFunction(const double t_position, const double) const
    {
        return t_position;
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularFunction
    //////////////////////////////////////////////////////////////////

    TabularFunction::TabularFunction(const std::vector<std::pair<double, double>> & values,
                                     const Variable property,
                                     const FenestrationCommon::Interpolator & interpolator) :
        IFunction(property),
        m_Curve(values),
        m_Interpolator(interpolator)
    {}

    TabularFunction::TabularFunction(const std::initializer_list<std::pair<double, double>> & list,
                                     const Variable property,
                                     const FenestrationCommon::Interpolator & interpolator) :
        IFunction(property),
        m_Curve(list),
        m_Interpolator(interpolator)
    {}


    double TabularFunction::evaluateFunction(const double t_position, const double) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](std::pair<double, double> val) {
            return val.first > t_position;
        });

        const auto points = getInterpolationPoints(it);

        return m_Interpolator.interpolate(points.first, points.second, t_position);
    }

    std::pair<std::pair<double, double>, std::pair<double, double>>
      TabularFunction::getInterpolationPoints(
        std::vector<std::pair<double, double>>::const_iterator & it) const
    {
        if(it == m_Curve.end())
        {
            --it;
        }
        const auto pt2 = *it;
        if(it != m_Curve.begin())
        {
            --it;
        }
        else
        {
            ++it;
        }

        const auto pt1 = *it;

        return std::make_pair(pt1, pt2);
    }

    double TabularFunction::maxX() const
    {
        return m_Curve.back().first;
    }

    double TabularFunction::maxY() const
    {
        return m_Curve.back().second;
    }

    double TabularFunction::minX() const
    {
        return m_Curve.front().first;
    }

    double TabularFunction::minY() const
    {
        return m_Curve.front().second;
    }

    const std::vector<std::pair<double, double>> & TabularFunction::getCurve() const
    {
        return m_Curve;
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivative
    //////////////////////////////////////////////////////////////////

    TabularDerivative::TabularDerivative(const std::vector<std::pair<double, double>> & values,
                                         const Variable property) :
        IFunction(property),
        m_Curve(values)
    {}

    TabularDerivative::TabularDerivative(
      const std::initializer_list<std::pair<double, double>> & list, const Variable property) :
        IFunction(property),
        m_Curve(list)
    {}

    double TabularDerivative::evaluateFunction(const double t_position, const double) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](std::pair<double, double> val) {
            return val.first > t_position;
        });
        const auto points = getInterpolationPoints(it);

        return (points.second.second - points.first.second)
               / (points.second.first - points.first.first);
    }

    std::pair<std::pair<double, double>, std::pair<double, double>>
      TabularDerivative::getInterpolationPoints(
        std::vector<std::pair<double, double>>::const_iterator & it) const
    {
        if(it == m_Curve.begin())
        {
            ++it;
        }
        auto pt2 = it == m_Curve.end() ? m_Curve.back() : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }
        if(*it == m_Curve.back())
        {
            --it;
        }
        auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

        return std::make_pair(pt1, pt2);
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivativeSmooth
    //////////////////////////////////////////////////////////////////

    TabularDerivativeSmooth::TabularDerivativeSmooth(
      const std::vector<std::pair<double, double>> & values, Variable property) :
        IFunction(property)
    {
        for(size_t i = 1u; i < values.size(); ++i)
        {
            const auto x1 = values[i - 1].first;
            const auto x2 = values[i].first;
            const auto y1 = values[i - 1].second;
            const auto y2 = values[i].second;
            const auto newX = (x1 + x2) / 2.0;
            const auto newY = (y2 - y1) / (x2 - x1);
            m_Curve.emplace_back(newX, newY);
        }
    }

    TabularDerivativeSmooth::TabularDerivativeSmooth(
      const std::initializer_list<std::pair<double, double>> & list, Variable property) :
        IFunction(property)
    {
        std::vector<std::pair<double, double>> helperVector{list};
        for(size_t i = 1u; i < helperVector.size(); ++i)
        {
            const auto x1 = helperVector[i - 1].first;
            const auto x2 = helperVector[i].first;
            const auto y1 = helperVector[i - 1].second;
            const auto y2 = helperVector[i].second;
            const auto newX = (x1 + x2) / 2.0;
            const auto newY = (y2 - y1) / (x2 - x1);
            m_Curve.emplace_back(newX, newY);
        }
    }

    double TabularDerivativeSmooth::evaluateFunction(const double t_position, const double) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](std::pair<double, double> val) {
            return val.first > t_position;
        });
        const auto points = getInterpolationPoints(it);

        const auto dy = points.second.second - points.first.second;
        const auto dx = points.second.first - points.first.first;

        auto result = points.first.second;
        if(dx != 0)
        {
            result += dy / dx * (t_position - points.first.first);
        }

        return result;
    }

    std::pair<std::pair<double, double>, std::pair<double, double>>
      TabularDerivativeSmooth::getInterpolationPoints(
        std::vector<std::pair<double, double>>::const_iterator & it) const
    {
        const auto pt2 = it == m_Curve.end() ? *std::prev(m_Curve.end()) : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }

        const auto pt1 = *it;

        return std::make_pair(pt1, pt2);
    }

    //////////////////////////////////////////////////////////////////
    ///  SuctionFunction
    //////////////////////////////////////////////////////////////////

    LiquidTransportationCurve::LiquidTransportationCurve(
      const std::vector<std::pair<double, double>> & values) :
        TabularFunction(values, Variable::water, FenestrationCommon::Interpolation::Logarithmic)
    {}

    LiquidTransportationCurve::LiquidTransportationCurve(
      const std::initializer_list<std::pair<double, double>> & list) :
        TabularFunction(list, Variable::water, FenestrationCommon::Interpolation::Logarithmic)
    {}

    std::pair<std::pair<double, double>, std::pair<double, double>>
      LiquidTransportationCurve::getInterpolationPoints(
        std::vector<std::pair<double, double>>::const_iterator & it) const
    {
        const auto pt2 = it == m_Curve.end() ? *std::prev(m_Curve.end()) : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }

        const auto pt1 = *it;

        return std::make_pair(pt1, pt2);
    }

    //////////////////////////////////////////////////////////////////
    ///  SaturationFunction
    //////////////////////////////////////////////////////////////////

    SaturationFunction::SaturationFunction() : IFunction(Variable::temperature)
    {}

    double SaturationFunction::evaluateFunction(const double t_position, const double) const
    {
        return saturationConcentrationAtTemperature(t_position);
    }

    //////////////////////////////////////////////////////////////////
    ///  Heat of evaporation
    //////////////////////////////////////////////////////////////////

    double HeatOfEvaporation::evaluateFunction(const double t_position, const double) const
    {
        return -(2500.8 - 2.36 * t_position + 0.016 * std::pow(t_position, 2)
                 - 0.00006 * std::pow(t_position, 3));
    }

    HeatOfEvaporation::HeatOfEvaporation() : IFunction(Variable::temperature)
    {}

    //////////////////////////////////////////////////////////////////
    ///  Phase change
    //////////////////////////////////////////////////////////////////
    PhaseChange::PhaseChange() : IFunction(Variable::temperature)
    {}

    double PhaseChange::evaluateFunction(const double t_position,
                                         const double t_PreviousTimestep) const
    {
        using Constants::FreezingPoint;
        using Constants::IcePoint;
        using Constants::EnthalpyOfFusion;

        auto result = 0.0;


        if(t_PreviousTimestep < IcePoint)
        {
            // Entire content of ice has melted.
            if(t_position >= FreezingPoint)
            {
                result = EnthalpyOfFusion;
            }

            // Only part of ice has melted and therefore linear interpolation is needed.
            if((t_position < FreezingPoint) && (t_position >= IcePoint))
            {
                result = EnthalpyOfFusion * (t_position - IcePoint) / (FreezingPoint - IcePoint);
            }
        }

        if(t_PreviousTimestep >= FreezingPoint)
        {
            // Entire content of water has froze.
            if(t_position < IcePoint)
            {
                result = -Constants::EnthalpyOfFusion;
            }

            // Only part of water froze and therefore linear interpolation is needed.
            if((t_position < FreezingPoint) && (t_position >= IcePoint))
            {
                result =
                  EnthalpyOfFusion * (t_position - FreezingPoint) / (FreezingPoint - IcePoint);
            }
        }

        return result;
    }

}   // namespace HygroThermFEM
