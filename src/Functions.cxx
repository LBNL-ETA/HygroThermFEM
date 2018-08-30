#include <algorithm>

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
        return evaluateFunction( state.getValue( m_Property ) );
    }

    //////////////////////////////////////////////////////////////////
    ///  IOperationFunction
    //////////////////////////////////////////////////////////////////

    double IOperation::value(const State & state) const
    {
        return m_Operator.at(m_Operation)(m_Function1->value(state), m_Function2->value(state));
    }

    IOperation::IOperation(std::shared_ptr<IValue> & t_Val1,
                           std::shared_ptr<IValue> & t_Val2,
                           Operation t_Operation) :
        m_Function1(t_Val1),
        m_Function2(t_Val2),
        m_Operation(t_Operation)
    {
        m_Operator[Operation::MULT] = [&](double a, double b) { return a * b; };
        m_Operator[Operation::DIV] = [&](double a, double b) { return a / b; };
        m_Operator[Operation::ADD] = [&](double a, double b) { return a + b; };
        m_Operator[Operation::SUB] = [&](double a, double b) { return a - b; };
    }

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    Constant::Constant(const double value) : IFunction(Property::temperature), m_Value(value)
    {}

    double Constant::evaluateFunction( const double ) const
    {
        return m_Value;
    }

    //////////////////////////////////////////////////////////////////
    ///  Operators
    //////////////////////////////////////////////////////////////////

    std::shared_ptr<MoisThermFEM::IValue> operator+(std::shared_ptr<IValue> & left,
                                                    std::shared_ptr<IValue> & right)
    {
        return std::make_shared<MoisThermFEM::IOperation>(
          left, right, MoisThermFEM::Operation::ADD);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator+(const double left,
                                                    std::shared_ptr<IValue> & right)
    {
        std::shared_ptr<IValue> aLeft(new MoisThermFEM::Constant(left));
        return std::make_shared<MoisThermFEM::IOperation>(
          aLeft, right, MoisThermFEM::Operation::ADD);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator+(std::shared_ptr<IValue> & left,
                                                    const double right)
    {
        return operator+(right, left);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator-(std::shared_ptr<IValue> & left,
                                                    std::shared_ptr<IValue> & right)
    {
        return std::make_shared<MoisThermFEM::IOperation>(
          left, right, MoisThermFEM::Operation::SUB);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator-(const double left,
                                                    std::shared_ptr<IValue> & right)
    {
        std::shared_ptr<MoisThermFEM::IValue> aLeft(new MoisThermFEM::Constant(left));
        return std::make_shared<MoisThermFEM::IOperation>(
          aLeft, right, MoisThermFEM::Operation::SUB);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator-(std::shared_ptr<IValue> & left,
                                                    const double right)
    {
        std::shared_ptr<MoisThermFEM::IValue> aRight(new MoisThermFEM::Constant(right));
        return std::make_shared<MoisThermFEM::IOperation>(
          left, aRight, MoisThermFEM::Operation::SUB);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator*(std::shared_ptr<IValue> & left,
                                                    std::shared_ptr<IValue> & right)
    {
        return std::make_shared<MoisThermFEM::IOperation>(
          left, right, MoisThermFEM::Operation::MULT);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator*(const double left,
                                                    std::shared_ptr<IValue> & right)
    {
        std::shared_ptr<MoisThermFEM::IValue> aLeft(new MoisThermFEM::Constant(left));
        return std::make_shared<MoisThermFEM::IOperation>(
          aLeft, right, MoisThermFEM::Operation::MULT);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator*(std::shared_ptr<IValue> & left,
                                                    const double right)
    {
        return operator*(right, left);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator/(std::shared_ptr<IValue> & left,
                                                    std::shared_ptr<IValue> & right)
    {
        return std::make_shared<MoisThermFEM::IOperation>(
          left, right, MoisThermFEM::Operation::DIV);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator/(const double left,
                                                    std::shared_ptr<IValue> & right)
    {
        std::shared_ptr<MoisThermFEM::IValue> aLeft(new MoisThermFEM::Constant(left));
        return std::make_shared<MoisThermFEM::IOperation>(
          aLeft, right, MoisThermFEM::Operation::DIV);
    }

    std::shared_ptr<MoisThermFEM::IValue> operator/(std::shared_ptr<IValue> & left,
                                                    const double right)
    {
        std::shared_ptr<MoisThermFEM::IValue> aRight(new MoisThermFEM::Constant(right));
        return std::make_shared<MoisThermFEM::IOperation>(
          left, aRight, MoisThermFEM::Operation::DIV);
    }

    //////////////////////////////////////////////////////////////////
    ///  Derivative
    //////////////////////////////////////////////////////////////////

    Derivative::Derivative(std::shared_ptr<IValue> & t_Value) : IValue(), m_Function(t_Value)
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

    //////////////////////////////////////////////////////////////////
    ///  TabularFunction
    //////////////////////////////////////////////////////////////////

    TabularFunction::TabularFunction(const std::vector<std::pair<double, double>> & values,
                                     Property property,
                                     FenestrationCommon::Interpolator interpolator) :
        IFunction(property),
        m_Curve(values),
        m_Interpolator(std::move(interpolator))
    {}

    TabularFunction::TabularFunction(std::initializer_list<std::pair<double, double>> & list,
                                     Property property,
                                     FenestrationCommon::Interpolator interpolator) :
        IFunction(property),
        m_Curve(list),
        m_Interpolator(std::move(interpolator))
    {}


    double TabularFunction::evaluateFunction( const double t_position ) const
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
                                         Property property) :
        IFunction(property),
        m_Curve(values)
    {}

    TabularDerivative::TabularDerivative(std::initializer_list<std::pair<double, double>> & list,
                                         Property property) :
        IFunction(property),
        m_Curve(list)
    {}

    double TabularDerivative::evaluateFunction( const double t_position ) const
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
                                     Property property,
                                     const FenestrationCommon::Interpolator & interpolator) :
        TabularFunction(values, property, interpolator)
    {}

    SuctionFunction::SuctionFunction(const std::initializer_list<std::pair<double, double>> & list,
                                     Property property,
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
        auto second = m_Curve.begin() + 1;
        auto pt2 = it == second ? m_Curve.front() : *it;
        if(it != m_Curve.begin())
        {
            --it;
        }

        auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

        return std::make_pair(pt1, pt2);
    }

    //////////////////////////////////////////////////////////////////
    ///  SaturationFunction
    //////////////////////////////////////////////////////////////////

    SaturationFunction::SaturationFunction(Property property, const double saturationCoefficient) :
        IFunction(property), m_SaturationCoefficient( saturationCoefficient )
    {}

    double SaturationFunction::evaluateFunction( const double t_position ) const
    {
        auto temp = 77.345 + 0.0057 * t_position - 7235.0 / t_position;
        temp = std::exp(temp);
        temp = temp / (461.4 * std::pow(t_position, m_SaturationCoefficient));
        return temp;
    }
}   // namespace MoisThermFEM
