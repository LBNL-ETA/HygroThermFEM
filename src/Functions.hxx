#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "Interpolator.hxx"
//#include "State.hxx"
//#include "Node2D.hxx"

/// Functions interface is used to build function that are used for matrix
/// building. Functions are stacked together to make full function that later
/// will be stored in FEM element.

namespace MoisThermFEM
{
    enum class Property;

    class State;
    class Node2D;
    class INodes;

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
        virtual double value(const Node2D & node) const = 0;
        virtual std::vector<double> values(const INodes & nodes) const;
        virtual ~IValue() = default;
    };

    using iValue = std::unique_ptr<IValue>;

    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    /// Interface for functions. Function must have defined property
    class IFunction : public IValue
    {
    public:
        IFunction(Property t_Property);

        double value(const Node2D & node) const override;

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
    template<class T, class U>
    class IOperation : public IValue
    {
    public:
        IOperation(const T t, const U s, const Operation & op) :
            m_Function1(std::move(t)),
            m_Function2(std::move(s)),
            m_Operation(op)
        {
            m_Operator[Operation::MULT] = [&](double a, double b) { return a * b; };
            m_Operator[Operation::DIV] = [&](double a, double b) { return a / b; };
            m_Operator[Operation::ADD] = [&](double a, double b) { return a + b; };
            m_Operator[Operation::SUB] = [&](double a, double b) { return a - b; };
        }

        double value(const Node2D & node) const override
        {
            return m_Operator.at(m_Operation)(m_Function1.value(node), m_Function2.value(node));
        };

    private:
        /// Functions can be shared between different operations and that is why it is
        /// necessary to share function
        const T m_Function1;
        const U m_Function2;

        const Operation m_Operation;

        /// This hold four basic operators (+, -. *. /) which is used to determine
        /// which function pointer is to be called
        std::map<Operation, std::function<double(double, double)>> m_Operator;
    };

    //////////////////////////////////////////////////////////////////
    ///  Operators for IValue
    //////////////////////////////////////////////////////////////////

    ////	operator+

    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator+(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::ADD);
    }

    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator+(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::ADD);
    }

    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator+(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::ADD);
    }

    ////	operator-

    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator-(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::SUB);
    }

    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator-(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::SUB);
    }

    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator-(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::SUB);
    }

    ////	operator*

    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator*(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::MULT);
    }

    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator*(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::MULT);
    }

    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator*(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::MULT);
    }

    ////	operator/

    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator/(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::DIV);
    }

    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator/(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::DIV);
    }

    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator/(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::DIV);
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
    ///  TabularFunction
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
    class SuctionCurve : public TabularFunction
    {
    public:
        SuctionCurve(const std::vector<std::pair<double, double>> & values,
                     const FenestrationCommon::Interpolator & interpolator =
                       FenestrationCommon::Interpolation::Logarithmic);

        SuctionCurve(const std::initializer_list<std::pair<double, double>> & list,
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
        SaturationFunction(const double saturationCoefficient = 9.2);

    private:
        double evaluateFunction(const double t_position) const override;
        const double m_SaturationCoefficient;
    };

}   // namespace MoisThermFEM
