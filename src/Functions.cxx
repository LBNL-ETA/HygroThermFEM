#include <algorithm>
#include <utility>
#include <vector>

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

    double heatOfEvaporation(const double temperature)
    {
        return (2500.8 - 2.36 * temperature + 0.016 * std::pow(temperature, 2)
                - 0.00006 * std::pow(temperature, 3))
               * 1000;
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

    double IFunction::value(const INode2D & node) const
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
    ///  TabularFunction1D
    //////////////////////////////////////////////////////////////////

    TabularFunction1D::TabularFunction1D(const std::vector<FenestrationCommon::point> & values,
                                         Variable property,
                                         FenestrationCommon::Interpolator interpolator) :
        IFunction(property),
        m_Curve(values),
        m_Interpolator(std::move(interpolator))
    {
        checkIfCurveIsSinglePoint();
    }

    TabularFunction1D::TabularFunction1D(
      const std::initializer_list<FenestrationCommon::point> & list,
      const Variable property,
      FenestrationCommon::Interpolator interpolator) :
        IFunction(property),
        m_Curve(list),
        m_Interpolator(std::move(interpolator))
    {
        checkIfCurveIsSinglePoint();
    }


    double TabularFunction1D::evaluateFunction(const double t_position, const double) const
    {
        auto it =
          std::find_if(m_Curve.begin(), m_Curve.end(), [&](const FenestrationCommon::point & val) {
              return val.x > t_position;
          });

        const auto points = getInterpolationPoints(it);

        return m_Interpolator.interpolate(points.first, points.second, t_position);
    }

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      TabularFunction1D::getInterpolationPoints(
        std::vector<FenestrationCommon::point>::const_iterator & it) const
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

    double TabularFunction1D::maxX() const
    {
        return m_Curve.back().x;
    }

    double TabularFunction1D::maxY() const
    {
        return m_Curve.back().y;
    }

    double TabularFunction1D::minX() const
    {
        return m_Curve.front().x;
    }

    double TabularFunction1D::minY() const
    {
        return m_Curve.front().y;
    }

    std::vector<FenestrationCommon::point> & TabularFunction1D::getCurve()
    {
        return m_Curve;
    }

    void TabularFunction1D::checkIfCurveIsSinglePoint()
    {
        if(m_Curve.size() == 1u)
        {
            auto point = m_Curve[0];
            point += FenestrationCommon::point(1, 0);
            m_Curve.emplace_back(point);
        }
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularFunction2D
    //////////////////////////////////////////////////////////////////

    TabularFunction2D::TabularFunction2D(
      const std::vector<FenestrationCommon::point> & firstValues,
      double firstTableMeasuredAt,
      Variable firstProperty,
      const std::vector<FenestrationCommon::point> & secondValues,
      double secondTableMeasureAt,
      Variable secondProperty,
      const FenestrationCommon::Interpolator & interpolator) :
        TabularFunction1D(secondValues, secondProperty, interpolator),
        m_FirstTable(firstValues, firstProperty, interpolator),
        m_FirstTableMeasuredAt{firstTableMeasuredAt},
        m_SecondTableMeasuredAt{secondTableMeasureAt}
    {
        const double tolerance{1e-6};
        const double y1{findValueAtPoint(firstValues, secondTableMeasureAt)};
        const double y2{findValueAtPoint(secondValues, firstTableMeasuredAt)};
        if(std::abs(y1 - y2) < tolerance)
        {
            m_CommonValueAtMeasuredTables = y1;
        }
        else
        {
            throw std::runtime_error("Values in two tables do not correspond to each other. Both "
                                     "tables must return same value at measured points.");
        }
    }

    double TabularFunction2D::value(const INode2D & node) const
    {
        const auto value1 = m_FirstTable.value(node);
        const auto value2 = TabularFunction1D::value(node);
        // This is simplified version of following result = m_CommonValueAtMeasuredTables + (value1
        // - m_CommonValueAtMeasuredTables) + (value2 - m_CommonValueAtMeasuredTables)
        return -m_CommonValueAtMeasuredTables + value1 + value2;
    }

    double TabularFunction2D::maxXFirstTable() const
    {
        return m_FirstTable.maxX();
    }

    double TabularFunction2D::maxYFirstTable() const
    {
        return m_FirstTable.maxY();
    }

    double TabularFunction2D::findValueAtPoint(const std::vector<FenestrationCommon::point> & table,
                                               const double value) const
    {
        auto it =
          std::find_if(table.begin(), table.end(), [&](const FenestrationCommon::point & val) {
              return val.x > value;
          });

        const auto points = getInterpolationPoints(it);

        return m_Interpolator.interpolate(points.first, points.second, value);
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivative
    //////////////////////////////////////////////////////////////////

    TabularDerivative::TabularDerivative(std::vector<FenestrationCommon::point> values,
                                         Variable property) :
        IFunction(property),
        m_Curve(std::move(values))
    {}

    TabularDerivative::TabularDerivative(
      const std::initializer_list<FenestrationCommon::point> & list, const Variable property) :
        IFunction(property),
        m_Curve(list)
    {}

    double TabularDerivative::evaluateFunction(const double t_position, const double) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](FenestrationCommon::point val) {
            return val.x > t_position;
        });
        const auto points = getInterpolationPoints(it);

        return (points.second.y - points.first.y) / (points.second.x - points.first.x);
    }

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      TabularDerivative::getInterpolationPoints(
        std::vector<FenestrationCommon::point>::const_iterator & it) const
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
      const std::vector<FenestrationCommon::point> & values, Variable property) :
        IFunction(property)
    {
        for(size_t i = 1u; i < values.size(); ++i)
        {
            const auto x1 = values[i - 1].x;
            const auto x2 = values[i].x;
            const auto y1 = values[i - 1].y;
            const auto y2 = values[i].y;
            const auto newX = (x1 + x2) / 2.0;
            const auto newY = (y2 - y1) / (x2 - x1);
            m_Curve.emplace_back(newX, newY);
        }
    }

    TabularDerivativeSmooth::TabularDerivativeSmooth(
      const std::initializer_list<FenestrationCommon::point> & list, Variable property) :
        IFunction(property)
    {
        std::vector<FenestrationCommon::point> helperVector{list};
        for(size_t i = 1u; i < helperVector.size(); ++i)
        {
            const auto x1 = helperVector[i - 1].x;
            const auto x2 = helperVector[i].x;
            const auto y1 = helperVector[i - 1].y;
            const auto y2 = helperVector[i].y;
            const auto newX = (x1 + x2) / 2.0;
            const auto newY = (y2 - y1) / (x2 - x1);
            m_Curve.emplace_back(newX, newY);
        }
    }

    double TabularDerivativeSmooth::evaluateFunction(const double t_position, const double) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](FenestrationCommon::point val) {
            return val.x > t_position;
        });
        const auto points = getInterpolationPoints(it);

        const auto dy = points.second.y - points.first.y;
        const auto dx = points.second.x - points.first.x;

        auto result = points.first.y;
        if(dx != 0)
        {
            result += dy / dx * (t_position - points.first.x);
        }

        return result;
    }

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      TabularDerivativeSmooth::getInterpolationPoints(
        std::vector<FenestrationCommon::point>::const_iterator & it) const
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
      const std::vector<FenestrationCommon::point> & vec) :
        TabularFunction1D(vec, Variable::water, FenestrationCommon::Interpolation::Logarithmic)
    {}

    LiquidTransportationCurve::LiquidTransportationCurve(
      const std::initializer_list<FenestrationCommon::point> & list) :
        TabularFunction1D(list, Variable::water, FenestrationCommon::Interpolation::Logarithmic)
    {}

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      LiquidTransportationCurve::getInterpolationPoints(
        std::vector<FenestrationCommon::point>::const_iterator & it) const
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
        return heatOfEvaporation(t_position);
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
