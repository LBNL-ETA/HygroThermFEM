#include <algorithm>
#include <cmath>

#include "FEMunique.hxx"
#include "Functions.hxx"
#include "State.hxx"
#include "Node2D.hxx"
#include "Common.hxx"

namespace MoisThermFEM
{
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
        return evaluateFunction(node.property(m_Property));
    }

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    Constant::Constant(const double value) : IFunction(Variable::temperature), m_Value(value)
    {}

    double Constant::evaluateFunction(const double) const
    {
        return m_Value;
    }

    //////////////////////////////////////////////////////////////////
    ///  State value
    //////////////////////////////////////////////////////////////////

    StateValue::StateValue(Variable property) : IFunction(property)
    {}

    double StateValue::evaluateFunction(const double t_position) const
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


    double TabularFunction::evaluateFunction(const double t_position) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](std::pair<double, double> val) {
            return val.first > t_position;
        });

        auto points = getInterpolationPoints(it);

        return m_Interpolator.interpolate(points.first, points.second, t_position);
    }

    std::pair<std::pair<double, double>, std::pair<double, double>>
      TabularFunction::getInterpolationPoints(
        std::vector<std::pair<double, double>>::const_iterator & it) const
    {
        auto pt2 = it == m_Curve.end() ? m_Curve.back() : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }
        auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

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

    double TabularDerivative::evaluateFunction(const double t_position) const
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
    ///  SuctionFunction
    //////////////////////////////////////////////////////////////////

    SuctionCurve::SuctionCurve(const std::vector<std::pair<double, double>> & values) :
        TabularFunction(values, Variable::humidity, FenestrationCommon::Interpolation::Logarithmic)
    {}

    SuctionCurve::SuctionCurve(const std::initializer_list<std::pair<double, double>> & list) :
        TabularFunction(list, Variable::humidity, FenestrationCommon::Interpolation::Logarithmic)
    {}

    std::pair<std::pair<double, double>, std::pair<double, double>>
      SuctionCurve::getInterpolationPoints(
        std::vector<std::pair<double, double>>::const_iterator & it) const
    {
        /// Suction curve takes care that first segment of curve always return value of first
        /// element.
        it == m_Curve.end() ? m_Curve.back() : *it;
        const auto second = std::next(m_Curve.begin());
        const auto pt2 = it == second ? m_Curve.front() : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }

        const auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

        return std::make_pair(pt1, pt2);
    }

    //////////////////////////////////////////////////////////////////
    ///  SaturationFunction
    //////////////////////////////////////////////////////////////////

    SaturationFunction::SaturationFunction(const double saturationCoefficient) :
        IFunction(Variable::temperature),
        m_SaturationCoefficient(saturationCoefficient)
    {}

    double SaturationFunction::evaluateFunction(const double t_position) const
    {
        const auto temperature = t_position + 273.15;
        auto temp = 77.345 + 0.0057 * temperature - 7235.0 / temperature;
        temp = std::exp(temp);
        temp = temp / (461.4 * std::pow(temperature, m_SaturationCoefficient));
        return temp;
    }

    //////////////////////////////////////////////////////////////////
    ///  Heat of evaporation
    //////////////////////////////////////////////////////////////////

    double HeatOfEvaporation::evaluateFunction(const double t_position) const
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

    double PhaseChange::evaluateFunction(double t_position) const
    {
        auto result = 0.0;
        if(t_position <= Constants::IcePoint)
        {
            result = Constants::EnthalpyOfFusion;
        }
        if((t_position < Constants::FreezingPoint) && (t_position >= Constants::IcePoint))
        {
            result = Constants::EnthalpyOfFusion * (t_position - Constants::IcePoint)
                /(Constants::FreezingPoint - Constants::IcePoint);
        }

        return result;
    }

}   // namespace MoisThermFEM
