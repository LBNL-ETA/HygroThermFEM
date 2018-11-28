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

    //! \brief Definition of interface used to perform various calculation(s) over state given in
    //! \brief node(s).
    //!
    //! When differential equation defines certain physical phenomenon contain functions that are
    //! dependent on some of state variables or they are simply constants, this class must be used
    //! as a parent.
    class IValue
    {
    public:
        //! Value at given node (redefined at every child class).
        virtual double value(const Node2D & node   //!< Node for which value is calculated for.
                             ) const = 0;

        //! Values at all nodes given with parameter INodes.
        virtual std::vector<double> values(const INodes & nodes   //!< Array of nodes
                                           ) const;
        virtual ~IValue() = default;
    };

    using iValue = std::unique_ptr<IValue>;

    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    //! \brief Interface for all functions that will be dependent on one of state variables.
    //!
    //! Any function that is dependent on any of state variables must be inherited from this class.
    //! Passed Property will be used to automatically read value of state variable and perform
    //! calculations.
    class IFunction : public IValue
    {
    public:
    	//! Basic constructor
        IFunction(
        	Property t_Property //!< Property for which function will be calculated.
        	);

    	//! Returns function evaluation for given node.
        double value(
        	const Node2D & node //!< Node at which function will be evaluated.
        	) const override;

    protected:
    	//! Interface for function definition. This is place where in inherited classes function
    	//! definitions will be stored.
        virtual double evaluateFunction(
        	const double t_position = 0 //!< Value at which function will be evaluated.
        		) const = 0;

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

	enum class Operation
	{
		MULT,
		DIV,
		ADD,
		SUB
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

    //////////////////////////////////////////////////////////////////
    ///  Heat of evaporation
    //////////////////////////////////////////////////////////////////

    class HeatOfEvaporation : public IFunction
    {
    public:
        HeatOfEvaporation();

    protected:
        double evaluateFunction(const double t_position) const override;
    };

}   // namespace MoisThermFEM
