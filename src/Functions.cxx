#include <algorithm>
#include <cmath>

#include "FEMunique.hxx"
#include "Functions.hxx"
#include "State.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM
{

	//////////////////////////////////////////////////////////////////
	///  IValue
	//////////////////////////////////////////////////////////////////

	std::vector< double > IValue::values(const INodes & nodes) const {
		std::vector< double > result(nodes.size(), 0);
		for(size_t i = 0u; i < nodes.size(); ++i)
		{
			result[i] = value(nodes[i].getState());
		}
		return result;
	}


    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    IFunction::IFunction(const Property t_Property) : m_Property(t_Property)
    {}

    double IFunction::value(const State & state) const
    {
        return evaluateFunction(state.getValue(m_Property));
    }

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    Constant::Constant(const double value) : IFunction(Property::temperature), m_Value(value)
    {}

    double Constant::evaluateFunction(const double) const
    {
        return m_Value;
    }

    //////////////////////////////////////////////////////////////////
    ///  State value
    //////////////////////////////////////////////////////////////////

    StateValue::StateValue(Property property) : IFunction(property)
    {}

    double StateValue::evaluateFunction(const double t_position) const
    {
        return t_position;
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularFunction
    //////////////////////////////////////////////////////////////////

    TabularFunction::TabularFunction(const std::vector<std::pair<double, double>> & values,
                                     const Property property,
                                     const FenestrationCommon::Interpolator & interpolator) :
        IFunction(property),
        m_Curve(values),
        m_Interpolator(interpolator)
    {}

    TabularFunction::TabularFunction(const std::initializer_list<std::pair<double, double>> & list,
                                     const Property property,
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

    double TabularFunction::max() const
    {
        return m_Curve.back().second;
    }

    double TabularFunction::min() const
    {
        return m_Curve.front().second;
    }

    std::vector<std::pair<double, double>> TabularFunction::getCurve() const
    {
        return m_Curve;
    }

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivative
    //////////////////////////////////////////////////////////////////

    TabularDerivative::TabularDerivative(const std::vector<std::pair<double, double>> & values,
                                         const Property property) :
        IFunction(property),
        m_Curve(values)
    {}

    TabularDerivative::TabularDerivative(
      const std::initializer_list<std::pair<double, double>> & list, const Property property) :
        IFunction(property),
        m_Curve(list)
    {}

    double TabularDerivative::evaluateFunction(const double t_position) const
    {
        auto it = std::find_if(m_Curve.begin(), m_Curve.end(), [&](std::pair<double, double> val) {
            return val.first > t_position;
        });
        auto points = getInterpolationPoints(it);

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

    SuctionFunction::SuctionFunction(const std::vector<std::pair<double, double>> & values,
                                     const Property property,
                                     const FenestrationCommon::Interpolator & interpolator) :
        TabularFunction(values, property, interpolator)
    {}

    SuctionFunction::SuctionFunction(const std::initializer_list<std::pair<double, double>> & list,
                                     const Property property,
                                     const FenestrationCommon::Interpolator & interpolator) :
        TabularFunction(list, property, interpolator)
    {}

    std::pair<std::pair<double, double>, std::pair<double, double>>
      SuctionFunction::getInterpolationPoints(
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

    SaturationFunction::SaturationFunction( const double saturationCoefficient ) :
        IFunction(Property::temperature),
        m_SaturationCoefficient(saturationCoefficient)
    {}

    double SaturationFunction::evaluateFunction(const double t_position) const
    {
        auto temp = 77.345 + 0.0057 * t_position - 7235.0 / t_position;
        temp = std::exp(temp);
        temp = temp / (461.4 * std::pow(t_position, m_SaturationCoefficient));
        return temp;
    }

}   // namespace MoisThermFEM
