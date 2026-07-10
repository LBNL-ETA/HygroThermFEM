#include <algorithm>
#include <utility>
#include <vector>
#include <stdexcept>

#include "Functions.hxx"
#include "Node2D.hxx"
#include "Material.hxx"
#include "Common.hxx"

namespace
{
    //! Finds the two bracketing points around iterator position for interpolation.
    //! Shared by TabularDerivativeSmooth and LiquidTransportationCurve.
    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      bracketingPoints(std::vector<FenestrationCommon::point>::const_iterator & iter,
                       const std::vector<FenestrationCommon::point> & table)
    {
        const auto pt2 = iter == table.end() ? *std::prev(table.end()) : *iter;
        if(iter != table.begin())
        {
            --iter;
        }
        const auto pt1 = *iter;
        return std::make_pair(pt1, pt2);
    }

    //! Computes midpoint derivative curve from a table of (x, y) points.
    //! Used by both TabularDerivativeSmooth constructors.
    std::vector<FenestrationCommon::point>
      computeDerivativeCurve(const std::vector<FenestrationCommon::point> & values)
    {
        std::vector<FenestrationCommon::point> curve;
        curve.reserve(values.size() - 1);
        for(size_t idx = 1u; idx < values.size(); ++idx)
        {
            const auto xLo = values[idx - 1].x;
            const auto xHi = values[idx].x;
            const auto yLo = values[idx - 1].y;
            const auto yHi = values[idx].y;
            curve.emplace_back((xLo + xHi) / 2.0, (yHi - yLo) / (xHi - xLo));
        }
        return curve;
    }
}   // anonymous namespace

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
        constexpr auto gasConstantForWaterVapor = 461.4;
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
        for(size_t idx = 0u; idx < nodes.size(); ++idx)
        {
            result[idx] = value(nodes[idx]);
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

        const auto points = getInterpolationPoints(it, m_Curve);

        return m_Interpolator.interpolate(points.first, points.second, t_position);
    }

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      TabularFunction1D::getInterpolationPoints(
          std::vector<FenestrationCommon::point>::const_iterator & it,
          const std::vector<FenestrationCommon::point> & table) const
    {
        if(it == table.end())
        {
            --it;
        }
        const auto pt2 = *it;
        if(table.size() != 1) {
            if (it != table.begin()) {
                --it;
            } else {
                ++it;
            }
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

    const std::vector<FenestrationCommon::point> & TabularFunction1D::getCurve() const
    {
        return m_Curve;
    }

    void TabularFunction1D::checkIfCurveIsSinglePoint()
    {
        if(m_Curve.size() == 1u)
        {
            m_Curve.emplace_back(m_Curve[0]);
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
        m_FirstTable(firstValues, firstProperty, interpolator)
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

        const auto points = getInterpolationPoints(it, table);

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
        const auto points = getInterpolationPoints(it, m_Curve);

        return (points.second.y - points.first.y) / (points.second.x - points.first.x);
    }

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      TabularDerivative::getInterpolationPoints(
          std::vector<FenestrationCommon::point>::const_iterator & it,
          const std::vector<FenestrationCommon::point> & table) const
    {
        if(it == table.begin())
        {
            ++it;
        }
        auto pt2 = it == table.end() ? table.back() : *it;
        if(it != table.begin())
        {
            --it;
        }
        if(*it == table.back())
        {
            --it;
        }
        auto pt1 = it == table.begin() ? table.front() : *it;

        return std::make_pair(pt1, pt2);
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivativeSmooth
    //////////////////////////////////////////////////////////////////

    TabularDerivativeSmooth::TabularDerivativeSmooth(
      const std::vector<FenestrationCommon::point> & values, Variable property) :
        IFunction(property),
        m_Curve(computeDerivativeCurve(values))
    {}

    TabularDerivativeSmooth::TabularDerivativeSmooth(
      const std::initializer_list<FenestrationCommon::point> & list, Variable property) :
        IFunction(property),
        m_Curve(computeDerivativeCurve(list))
    {}

    double TabularDerivativeSmooth::evaluateFunction(const double t_position, const double) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](FenestrationCommon::point val) {
            return val.x > t_position;
        });
        const auto points = getInterpolationPoints(it, m_Curve);

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
          std::vector<FenestrationCommon::point>::const_iterator & iter,
          const std::vector<FenestrationCommon::point> & table) const
    {
        return bracketingPoints(iter, table);
    }

    //////////////////////////////////////////////////////////////////
    ///  SorptionSecantCapacity
    //////////////////////////////////////////////////////////////////

    SorptionSecantCapacity::SorptionSecantCapacity(
      const std::vector<FenestrationCommon::point> & sorptionCurve, const Variable property) :
        TabularFunction1D(sorptionCurve, property)
    {}

    double SorptionSecantCapacity::evaluateFunction(const double t_position,
                                                    const double t_previousTimestep) const
    {
        // w(phi): the sorption curve, via the base class' linear interpolation.
        const auto w = [this](const double phi) {
            return TabularFunction1D::evaluateFunction(phi, phi);
        };

        const double dphi = t_position - t_previousTimestep;
        constexpr double minDelta = 1e-12;
        if(std::abs(dphi) > minDelta)
        {
            // Secant slope over the step -> mass-conservative capacity.
            return (w(t_position) - w(t_previousTimestep)) / dphi;
        }

        // Secant undefined (no change since previous timestep): finite-difference tangent.
        constexpr double h = 1e-6;
        return (w(t_position + h) - w(t_position - h)) / (2.0 * h);
    }

    //////////////////////////////////////////////////////////////////
    ///  SuctionFunction
    //////////////////////////////////////////////////////////////////

    LiquidTransportationCurve::LiquidTransportationCurve(
      const std::vector<FenestrationCommon::point> & vec, const IMaterial & material) :
        TabularFunction1D(vec, Variable::water, FenestrationCommon::Interpolation::Logarithmic),
        m_Material(material)
    {}

    LiquidTransportationCurve::LiquidTransportationCurve(
      const std::initializer_list<FenestrationCommon::point> & list, const IMaterial & material) :
        TabularFunction1D(list, Variable::water, FenestrationCommon::Interpolation::Logarithmic),
        m_Material(material)
    {}

    double LiquidTransportationCurve::value(const INode2D & node) const
    {
        // Resolve `w` through the bound material instead of through the node's
        // averaged `Variable::water` property. At a material interface the node's
        // averaged water content is not consistent with this material's curve at
        // the node's humidity, which produces a non-uniform per-element coefficient
        // and breaks the property `K * U_uniform = 0` of the assembled stiffness.
        // See class comment in Functions.hxx for the full reasoning.
        const double materialWater =
          m_Material.waterContent(node).content(WaterContent::Water);
        return evaluateFunction(materialWater, materialWater);
    }

    std::pair<FenestrationCommon::point, FenestrationCommon::point>
      LiquidTransportationCurve::getInterpolationPoints(
          std::vector<FenestrationCommon::point>::const_iterator & iter,
          const std::vector<FenestrationCommon::point> & table) const
    {
        return bracketingPoints(iter, table);
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
