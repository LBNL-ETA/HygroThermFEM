#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "Interpolator.hxx"
#include "VectorOperators.hxx"
#include "State.hxx"

/// Functions interface is used to build function that are used for matrix
/// building. Functions are stacked together to make full function that later
/// will be stored in FEM element.

namespace MoisThermFEM
{
    enum class Property;

    class State;

    enum class Operation
    {
        MULT,
        DIV,
        ADD,
        SUB
    };

    class IValue
    {
    public:
        virtual double value(const State & state) const = 0;
        virtual ~IValue() = default;
    };

    using iValue = std::unique_ptr<IValue>;

    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    /// Interface for functions
    class IFunction : public IValue
    {
    public:
        IFunction(Property t_Property);

        virtual double value(const State & state) const override;

    protected:
        virtual double evaluateFunction(const double t_position = 0) const = 0;

        /// Property that is used to calculate function value. It is extracted from current
        /// domain (material) point.
        const Property m_Property;
    };

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    /// Simple constant curve.
    class Constant : public IFunction
    {
    public:
        Constant(const double value);

    private:
        double evaluateFunction(const double t_position) const override;

        double m_Value;
    };

    //////////////////////////////////////////////////////////////////
    ///  IOperation
    //////////////////////////////////////////////////////////////////

    /// Entire class is used to mimic operator functions so that FEM functions can be
    /// written as ordinary equations.
    template<class T, class U, class OPERATION>
    class IOperation : public IValue
    {
    public:
        IOperation(const T & t, const U & s, const OPERATION & op) :
            m_Function1(t),
            m_Function2(s),
            m_Operation(op)
        {
            m_Operator[Operation::MULT] = [&](double a, double b) { return a * b; };
            m_Operator[Operation::DIV] = [&](double a, double b) { return a / b; };
            m_Operator[Operation::ADD] = [&](double a, double b) { return a + b; };
            m_Operator[Operation::SUB] = [&](double a, double b) { return a - b; };
        }

        double value(const State & state) const override
        {
            return m_Operator.at(m_Operation)(m_Function1.value(state), m_Function2.value(state));
        };

    private:
        /// Functions can be shared between different operations and that is why it is
        /// necessary to share function
        const T m_Function1;
        const U m_Function2;

        const OPERATION m_Operation;

        /// This hold four basic operators (+, -. *. /) which is used to determine
        /// which function pointer is to be called
        std::map<Operation, std::function<double(double, double)>> m_Operator;
    };

    //////////////////////////////////////////////////////////////////
    ///  Operators
    //////////////////////////////////////////////////////////////////

    template<class T>
    IOperation<T, Constant, Operation> operator+(const T & t, double u)
    {
        Constant con{u};
        return IOperation<T, Constant, Operation>(t, con, Operation::ADD);
    }

    template<class T>
    IOperation<Constant, T, Operation> operator+(double u, const T & t)
    {
        Constant con{u};
        return IOperation<Constant, T, Operation>(con, t, Operation::ADD);
    }

    template<class T, class U>
    IOperation<T, U, Operation> operator+(T const & t, U const & u)
    {
        return IOperation<T, U, Operation>(t, u, Operation::ADD);
    }

    template<class T>
    IOperation<T, Constant, Operation> operator-(const T & t, double u)
    {
        Constant con{u};
        return IOperation<T, Constant, Operation>(t, con, Operation::SUB);
    }

    template<class T>
    IOperation<Constant, T, Operation> operator-(double u, const T & t)
    {
        Constant con{u};
        return IOperation<Constant, T, Operation>(con, t, Operation::SUB);
    }

    template<class T, class U>
    IOperation<T, U, Operation> operator-(T const & t, U const & u)
    {
        return IOperation<T, U, Operation>(t, u, Operation::SUB);
    }

    template<class T>
    IOperation<T, Constant, Operation> operator*(const T & t, double u)
    {
        Constant con{u};
        return IOperation<T, Constant, Operation>(t, con, Operation::MULT);
    }

    template<class T>
    IOperation<Constant, T, Operation> operator*(double u, const T & t)
    {
        Constant con{u};
        return IOperation<Constant, T, Operation>(con, t, Operation::MULT);
    }

    template<class T, class U>
    IOperation<T, U, Operation> operator*(T const & t, U const & u)
    {
        return IOperation<T, U, Operation>(t, u, Operation::MULT);
    }

    template<class T>
    IOperation<T, Constant, Operation> operator/(const T & t, double u)
    {
        Constant con{u};
        return IOperation<T, Constant, Operation>(t, con, Operation::DIV);
    }

    template<class T>
    IOperation<Constant, T, Operation> operator/(double u, const T & t)
    {
        Constant con{u};
        return IOperation<Constant, T, Operation>(con, t, Operation::DIV);
    }

    template<class T, class U>
    IOperation<T, U, Operation> operator/(T const & t, U const & u)
    {
        return IOperation<T, U, Operation>(t, u, Operation::DIV);
    }

    //////////////////////////////////////////////////////////////////
    ///  State value
    //////////////////////////////////////////////////////////////////

    class StateValue : public IFunction
    {
    public:
        StateValue(Property property);

    private:
        double evaluateFunction(const double t_position) const override;
    };

    //////////////////////////////////////////////////////////////////
    ///  Derivative
    //////////////////////////////////////////////////////////////////
    template<class T>
    class Derivative : public IValue
    {
    public:
        Derivative(const T & t_Value) : IValue(), m_Function(t_Value)
        {}

        double value(const State & state) const override
        {
            double val1 = m_Function.value(state);

            // small depends on exact number that we are calculating
            const double small = val1 != 0 ? val1 / 1e5 : 1e-5;
            const State smallIncrease(state + State(small, small, small, 0));

            double val2 = m_Function.value(smallIncrease);
            return (val2 - val1) / small;
        }

    private:
        const T m_Function;
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularFunctions
    //////////////////////////////////////////////////////////////////

    /// Interface for classic tabular curve. There are different interpolation
    /// strategies and this is base class for all of them.
    class TabularFunction : public IFunction
    {
    public:
        TabularFunction(const std::vector<std::pair<double, double>> & values,
                        const Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Linear);

        TabularFunction(const std::initializer_list<std::pair<double, double>> & list,
                        const Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Linear);

        double max() const;

        double min() const;

        std::vector<std::pair<double, double>> getCurve() const;

    protected:
        std::vector<std::pair<double, double>> m_Curve;
        FenestrationCommon::Interpolator m_Interpolator;

        double evaluateFunction(const double t_position) const override;

        virtual std::pair<std::pair<double, double>, std::pair<double, double>>
          getInterpolationPoints(std::vector<std::pair<double, double>>::const_iterator & it) const;
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivative
    //////////////////////////////////////////////////////////////////

    /// This class is different from ordinary derivative because it extends over the
    /// limits. This is important in iterations when first derivative really needs
    /// to be evaluated outside of limits or convergence will produce incorrect
    /// results (sorption curve is good example).
    class TabularDerivative : public IFunction
    {
    public:
        TabularDerivative(const std::vector<std::pair<double, double>> & values, Property property);

        TabularDerivative(const std::initializer_list<std::pair<double, double>> & list,
                          Property property);

    protected:
        std::vector<std::pair<double, double>> m_Curve;

        double evaluateFunction(const double t_position) const override;

        virtual std::pair<std::pair<double, double>, std::pair<double, double>>
          getInterpolationPoints(std::vector<std::pair<double, double>>::const_iterator & it) const;
    };

    //////////////////////////////////////////////////////////////////
    ///  SuctionFunction
    //////////////////////////////////////////////////////////////////

    /// Class that behaves like suction curve. It is standard (linear or
    /// logarithmic) interpolation except for the results in first range where curve
    /// will return constant value equal to the first point
    class SuctionFunction : public TabularFunction
    {
    public:
        SuctionFunction(const std::vector<std::pair<double, double>> & values,
                        Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Logarithmic);

        SuctionFunction(const std::initializer_list<std::pair<double, double>> & list,
                        Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Logarithmic);

    protected:
        std::pair<std::pair<double, double>, std::pair<double, double>> getInterpolationPoints(
          std::vector<std::pair<double, double>>::const_iterator & it) const override;
    };

    //////////////////////////////////////////////////////////////////
    ///  SaturationFunction
    //////////////////////////////////////////////////////////////////

    /// Simple constant curve.
    class SaturationFunction : public IFunction
    {
    public:
        SaturationFunction(Property property, double saturationCoefficient = 9.2);

    private:
        double evaluateFunction(const double t_position) const override;
        const double m_SaturationCoefficient;
    };

}   // namespace MoisThermFEM
