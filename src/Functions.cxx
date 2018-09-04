#include <algorithm>
#include <cmath>

#include "FEMunique.hxx"
#include "Functions.hxx"
#include "State.hxx"

namespace MoisThermFEM
{
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
    ///  IOperationFunction
    //////////////////////////////////////////////////////////////////

    double IOperation::value(const State & state) const
    {
        return m_Operator.at(m_Operation)(m_Function1->value(state), m_Function2->value(state));
    }

    IOperation::IOperation(iValue t_Val1, iValue t_Val2, const Operation t_Operation) :
        m_Function1(t_Val1->clone()),
        m_Function2(t_Val2->clone()),
        m_Operation(t_Operation)
    {
        m_Operator[Operation::MULT] = [&](double a, double b) { return a * b; };
        m_Operator[Operation::DIV] = [&](double a, double b) { return a / b; };
        m_Operator[Operation::ADD] = [&](double a, double b) { return a + b; };
        m_Operator[Operation::SUB] = [&](double a, double b) { return a - b; };
    }

    iValue IOperation::clone() const
    {
        return fem::make_unique<IOperation>(
          this->m_Function1->clone(), this->m_Function2->clone(), m_Operation);
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

    std::unique_ptr<Constant> Constant::create(const double value)
    {
        return std::unique_ptr<Constant>(new Constant(value));
    }

    iValue Constant::clone() const
    {
        return fem::make_unique<Constant>(*this);
    }

    //////////////////////////////////////////////////////////////////
    ///  Operators
    //////////////////////////////////////////////////////////////////

    iValue operator+(iValue & left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(left->clone(), right->clone(), Operation::ADD));
    }

    iValue operator+(const double left, iValue & right)
    {
        iValue aLeft = Constant::create(left);
        return std::unique_ptr<IOperation>(
          new IOperation(aLeft->clone(), right->clone(), Operation::ADD));
    }

    iValue operator+(iValue & left, const double right)
    {
        return operator+(right, left);
    }

    iValue operator-(iValue & left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(left->clone(), right->clone(), Operation::SUB));
    }

    iValue operator-(const double left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(Constant::create(left), right->clone(), Operation::SUB));
    }

    iValue operator-(iValue & left, const double right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(left->clone(), Constant::create(right), Operation::SUB));
    }

    iValue operator*(iValue & left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(left->clone(), right->clone(), Operation::MULT));
    }

    iValue operator*(const double left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(Constant::create(left), right->clone(), MoisThermFEM::Operation::MULT));
    }

    iValue operator*(iValue & left, const double right)
    {
        return operator*(right, left);
    }

    iValue operator/(iValue & left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(left->clone(), right->clone(), Operation::DIV));
    }

    iValue operator/(const double left, iValue & right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(Constant::create(left), right->clone(), Operation::DIV));
    }

    iValue operator/(iValue & left, const double right)
    {
        return std::unique_ptr<IOperation>(
          new IOperation(left->clone(), Constant::create(right), Operation::DIV));
    }

    //////////////////////////////////////////////////////////////////
    ///  Derivative
    //////////////////////////////////////////////////////////////////

    Derivative::Derivative(const iValue & t_Value) : IValue(), m_Function(t_Value->clone())
    {}

    double Derivative::value(const State & state) const
    {
        double val1 = m_Function->value(state);

        // small depends on exact number that we are calculating
        const double small = val1 != 0 ? val1 / 1e5 : 1e-5;
        const State smallIncrease(state + State(small, small, small));

        double val2 = m_Function->value(smallIncrease);
        return (val2 - val1) / small;
    }

    std::unique_ptr<Derivative> Derivative::create(const iValue & t_Function)
    {
        return std::unique_ptr<Derivative>(new Derivative(t_Function));
    }

    iValue Derivative::clone() const
    {
        return create(m_Function->clone());
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

    std::unique_ptr<TabularFunction>
      TabularFunction::create(const std::initializer_list<std::pair<double, double>> & list,
                              const Property property,
                              const FenestrationCommon::Interpolator & interpolator)
    {
        return std::unique_ptr<TabularFunction>(new TabularFunction(list, property, interpolator));
    }

    std::unique_ptr<TabularFunction>
      TabularFunction::create(const std::vector<std::pair<double, double>> & values,
                              const Property property,
                              const FenestrationCommon::Interpolator & interpolator)
    {
        return std::unique_ptr<TabularFunction>(
          new TabularFunction(values, property, interpolator));
    }

    iValue TabularFunction::clone() const
    {
        return fem::make_unique<TabularFunction>(*this);
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
		const std::initializer_list< std::pair< double, double>> & list,
		const Property property ) :
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

    std::unique_ptr<TabularDerivative>
      TabularDerivative::create(const std::vector<std::pair<double, double>> & values,
                                const Property property)
    {
        return std::unique_ptr<TabularDerivative>(new TabularDerivative(values, property));
    }

    std::unique_ptr<TabularDerivative>
      TabularDerivative::create( const std::initializer_list< std::pair< double, double>> & list,
								 Property property )
    {
        return std::unique_ptr<TabularDerivative>(new TabularDerivative(list, property));
    }

    iValue TabularDerivative::clone() const
    {
        return fem::make_unique<TabularDerivative>(*this);
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
        const auto second = m_Curve.begin() + 1;
        const auto pt2 = it == second ? m_Curve.front() : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }

        const auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

        return std::make_pair(pt1, pt2);
    }

    std::unique_ptr<SuctionFunction>
      SuctionFunction::create(const std::vector<std::pair<double, double>> & values,
                              const Property property,
                              const FenestrationCommon::Interpolator & interpolator)
    {
        return std::unique_ptr<SuctionFunction>(
          new SuctionFunction(values, property, interpolator));
    }

    std::unique_ptr<SuctionFunction>
      SuctionFunction::create(const std::initializer_list<std::pair<double, double>> & list,
                              const Property property,
                              const FenestrationCommon::Interpolator & interpolator)
    {
        return std::unique_ptr<SuctionFunction>(new SuctionFunction(list, property, interpolator));
    }

    iValue SuctionFunction::clone() const
    {
        return fem::make_unique<SuctionFunction>(*this);
    }

    //////////////////////////////////////////////////////////////////
    ///  SaturationFunction
    //////////////////////////////////////////////////////////////////

    SaturationFunction::SaturationFunction(const Property property, const double saturationCoefficient) :
        IFunction(property),
        m_SaturationCoefficient(saturationCoefficient)
    {}

    double SaturationFunction::evaluateFunction(const double t_position) const
    {
        auto temp = 77.345 + 0.0057 * t_position - 7235.0 / t_position;
        temp = std::exp(temp);
        temp = temp / (461.4 * std::pow(t_position, m_SaturationCoefficient));
        return temp;
    }

    std::unique_ptr<SaturationFunction> SaturationFunction::create(const Property property,
                                                                   double saturationCoefficient)
    {
        return std::unique_ptr<SaturationFunction>(
          new SaturationFunction(property, saturationCoefficient));
    }

    iValue SaturationFunction::clone() const
    {
        return fem::make_unique<SaturationFunction>(*this);
    }
}   // namespace MoisThermFEM
