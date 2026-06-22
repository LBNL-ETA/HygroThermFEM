#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "Interpolator.hxx"
#include "Point.hxx"

/// Functions interface is used to build function that are used for matrix
/// building. Functions are stacked together to make full function that later
/// will be stored in FEM element.

namespace HygroThermFEM
{
    enum class Variable;

    struct State;

    class INode2D;

    class INodes;

    class IMaterial;

    //! \brief Calculates vapor pressure at given temperature
    //!
    //! \param temperature Temperature in Celsius
    //! \return vapor pressure [Pa]
    double vaporPressureAtTemperature(double temperature);

    //! \brief Saturation function
    //!
    //! \param temperature Temperature in Celsius
    //! \return water content [kg/m^3]
    double saturationConcentrationAtTemperature(double temperature);

    //! Heat of evaporation
    //!
    //! \param temperature Temperature in Celsius
    //! \return Heat of evaporation
    double heatOfEvaporation(double temperature);

    //! \brief Definition of interface used to perform various calculation(s) over state given in
    //! \brief node(s).
    //!
    //! When differential equation defines certain physical phenomenon contain functions that are
    //! dependent on some of state variables or they are simply constants, this class must be used
    //! as a parent.
    class IValue
    {
    public:
        //! \brief Value at given node (redefined at every child class).
        //!
        //! \param node Node for which value is calculated for.
        //! \return State variable value for given node
        virtual double value(const INode2D & node) const = 0;

        //! \brief Values at all nodes given with parameter INodes.
        //!
        //! \param nodes Array of nodes for which values are being requested
        //! \return Array of state variable values for requested nodes
        virtual std::vector<double> values(const INodes & nodes) const;

        //! \brief Missing default destructor in abstract class can cause memory leaks.
        virtual ~IValue() = default;
    };

    using iValue = std::unique_ptr<IValue>;

    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    //! \brief Interface for all functions that will be dependent on one of state variables.
    //!
    //! Any function that is dependent on any of state variables must be inherited from this class.
    //! Passed Variable will be used to automatically read value of state variable and perform
    //! calculations.
    class IFunction : public IValue
    {
    public:
        //! \brief Basic constructor
        //!
        //! \param t_Property Variable for which function will be calculated.
        IFunction(Variable t_Property);

        //! \brief Returns function evaluation for given node.
        //!
        //! \param node
        //! \return Node at which function will be evaluated.
        virtual double value(const INode2D & node) const override;

    protected:
        //! \brief Interface for function definition. This is place where in inherited classes
        //! function \brief definitions will be stored.
        //!
        //! \param t_position Value at which function will be evaluated.
        //! \param t_previousTimestep Value from previous timestep.
        //! \return evaluated state variable value.
        virtual double evaluateFunction(double t_position = 0,
                                        double t_previousTimestep = 0) const = 0;

        /// Variable that is used to calculate function value. It is extracted from current
        /// domain (material) point.
        const Variable m_Property;
    };

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    //! \brief Simple constant variable that will be used in function definitions for finite
    //! elements.
    class Constant : public IFunction
    {
    public:
        Constant(double value);

    private:
        double evaluateFunction(double t_position, double t_previousTimestep) const override;

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

    //! \brief It is used to support operation with functions.
    //!
    //! Functions that are child of IValue are used directly in differential equations.
    //! Those functions can be stacked with ordinary operations. This class is used to support
    //! those operations.
    template<class T, class U>
    class IOperation : public IValue
    {
    public:
        //! Constructor that accept two operands and operation.
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

        //! Returns value of operation.
        double value(const INode2D & node) const override
        {
            return m_Operator.at(m_Operation)(m_Function1.value(node), m_Function2.value(node));
        }

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

    //! template operator + for any class derived from IValue. It accepts two operands.
    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator+(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::ADD);
    }

    //! template operator + for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator+(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::ADD);
    }

    //! template operator + for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator+(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::ADD);
    }

    ////	operator-

    //! template operator - for any class derived from IValue. It accepts two operands.
    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator-(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::SUB);
    }

    //! template operator - for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator-(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::SUB);
    }

    //! template operator - for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator-(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::SUB);
    }

    ////	operator*

    //! template operator * for any class derived from IValue. It accepts two operands.
    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator*(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::MULT);
    }

    //! template operator * for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator*(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::MULT);
    }

    //! template operator * for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename U>
    typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
      operator*(const double & t, const U & u)
    {
        Constant con{t};
        return IOperation<Constant, U>(con, u, Operation::MULT);
    }

    ////	operator/

    //! template operator / for any class derived from IValue. It accepts two operands.
    template<typename T, typename U>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
      operator/(const T & t, const U & u)
    {
        return IOperation<T, U>(t, u, Operation::DIV);
    }

    //! template operator / for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
    template<typename T>
    typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
      operator/(const T & t, const double & u)
    {
        Constant con{u};
        return IOperation<T, Constant>(t, con, Operation::DIV);
    }

    //! template operator / for any class derived from IValue. It accepts two operands one of which
    //! is type of double.
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

    //! \brief Class that simply returns value of state variable.
    //!
    //! In some of differential equations, value of current state variable can depend on value of
    //! state variable previously calculated. For example, temperature distribution will be
    //! dependent on water content.
    class StateValue : public IFunction
    {
    public:
        //! \brief Constructor
        //!
        //! \param property Variable that StateValue represents.
        StateValue(Variable property);

    private:
        //! \brief Inherited function evaluation for current property
        //!
        //! \param t_position Position for which state value will be evaluated at. In this case it
        //! is simply state variable value itself.
        //! \param t_previousTimestep Value at previous timestep (Not used in this case).
        //! \return Value at given postion (Equal to position in this case.
        double evaluateFunction(double t_position, double t_previousTimestep) const override;
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularFunction1D
    //////////////////////////////////////////////////////////////////

    //! \brief Interface for classic tabular curve.
    //!
    //! Tabular function is made to simply return adequate value from table. Search function is
    //! performed over first column and return value is calculated from second column.
    class TabularFunction1D : public IFunction
    {
    public:
        //! \brief Table construction from standard vector.
        //!
        //! \param values Vector of x,y pairs that will go into table.
        //! \param property State variable that represent value type in first column.
        //! \param interpolator Interpolation strategy used for in-between values.
        TabularFunction1D(const std::vector<FenestrationCommon::point> &values,
                          Variable property,
                          FenestrationCommon::Interpolator interpolator =
                          FenestrationCommon::Interpolation::Linear);

        //! \brief Table construction from initializer list.
        //!
        //! \param list Initializer list of x,y pairs that will go into table.
        //! \param property Variable FenestrationCommon::Interpolator n.
        //! \param interpolator Interpolation strategy used for in-between values.
        TabularFunction1D(const std::initializer_list<FenestrationCommon::point> & list,
                          Variable property,
                          FenestrationCommon::Interpolator interpolator =
                            FenestrationCommon::Interpolation::Linear);

        //! \brief Returns maximum value for first column.
        double maxX() const;

        //! \brief Returns maximum value from second column.
        double maxY() const;

        //! \brief FenestrationCommon::Interpolator
        double minX() const;

        //! \brief Returns minimum value from second column.
        double minY() const;

        //! \brief Returns tabular values as standard vector of pairs.
        std::vector<FenestrationCommon::point> & getCurve();

        //! \brief Returns tabular values as standard vector of pairs (read-only).
        const std::vector<FenestrationCommon::point> & getCurve() const;

    protected:
        std::vector<FenestrationCommon::point> m_Curve;
        FenestrationCommon::Interpolator m_Interpolator;

        //! \brief Function to calculate value from given table.
        virtual double evaluateFunction(double t_position,
                                        double t_previousTimestep) const override;

        //! \brief Helper function that returns two closest points for interpolation.
        virtual std::pair<FenestrationCommon::point, FenestrationCommon::point>
          getInterpolationPoints(std::vector<FenestrationCommon::point>::const_iterator & it, const std::vector<FenestrationCommon::
                                 point> & table) const;

        //! \brief Expands table into two points if only one point is provided.
        void checkIfCurveIsSinglePoint();
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularFunction2D
    //////////////////////////////////////////////////////////////////

    //! \brief Interface for tabular function in 2D space.
    //!
    //! Some properties are provided with two dimensional table. Good example is thermal
    //! conductivity that is both temperature and moisture dependent.
    class TabularFunction2D : public TabularFunction1D
    {
    public:
        //! \brief Constructor of 2D table
        //! 2D table will be constructed of two tables that will form two dimensional table.
        //!
        //! \param firstValues Values for first table dependency
        //! \param firstTableMeasuredAt At which value is first table measured at. If for example, first table is
        //! temperature dependent and second table is humidity dependent, then this value should represent.
        //! at which humidity first table values are measured at.
        //! \param firstProperty Property of first table.
        //! \param secondValues Values for second table dependency.
        //! \param secondTableMeasureAt Value at which second table is measure at.
        //! \param secondProperty Property of second table.
        //! \param interpolator Interpolation strategy.
        TabularFunction2D(const std::vector<FenestrationCommon::point> &firstValues,
                          double firstTableMeasuredAt,
                          Variable firstProperty,
                          const std::vector<FenestrationCommon::point> &secondValues,
                          double secondTableMeasureAt,
                          Variable secondProperty,
                          const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Linear);

        double value(const INode2D & node) const override;

        double maxXFirstTable() const;
        double maxYFirstTable() const;

    private:

        //! \brief Function to find y value from vector of points
        //!
        //! \param table Table of x-y points
        //! \param value x-value in table for which y-value needs to be evaluated.
        //! \return y-value for corresponding passed x-value
        double findValueAtPoint(const std::vector<FenestrationCommon::point> & table, double value) const;

        TabularFunction1D m_FirstTable;

        // Both tables need to match resulting value at measured properties
        double m_CommonValueAtMeasuredTables;
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivative
    //////////////////////////////////////////////////////////////////

    //! \brief Estimates tabular derivative.
    //!
    //! This class is different from ordinary derivative because it extends over the
    //! limits. This is important in iterations when first derivative really needs
    //! to be evaluated outside of limits or convergence will produce incorrect
    //! results (sorption curve is good example).
    class TabularDerivative : public IFunction
    {
    public:
        //! \brief Construction of tabular derivative from standard vector values.
        //!
        //! \param values Points that construct the table.
        //! \param property Variable that represent value type in x coordinate.
        TabularDerivative(std::vector<FenestrationCommon::point> values, Variable property);

        //! \brief Construction of tabular derivative from initializer_list.
        //!
        //! \param list Initializer list used to construct tabular derivative.
        //! \param property Variable that represent value type in first column.
        TabularDerivative(const std::initializer_list<FenestrationCommon::point> & list,
                          Variable property);

    protected:
        std::vector<FenestrationCommon::point> m_Curve;

        //! \brief Evaluation function for given class (override from base class).
        //!
        //! \param t_position value for which function needs to be evaluated.
        //! \param t_previousTimestep Value in previous timestep.
        //! \return Value for given function at requested position.
        double evaluateFunction(double t_position, double t_previousTimestep) const override;

        //! \brief Helper function that returns two closest points for interpolation.
        virtual std::pair<FenestrationCommon::point, FenestrationCommon::point>
          getInterpolationPoints(std::vector<FenestrationCommon::point>::const_iterator & it,
                                 const std::vector<FenestrationCommon::point> & table) const;
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularDerivativeSmooth
    //////////////////////////////////////////////////////////////////

    //! \brief Estimates tabular derivative with values being smoothed between different values.
    //!
    //! When providing table with x-y values, calculating simple derivative will cause step function
    //! between certain values. This class will provide smooth transition between different values.
    class TabularDerivativeSmooth : public IFunction
    {
    public:
        //! \brief Construction of tabular derivative from standard vector values.
        //!
        //! \param values Points that construct the table
        //! \param property Variable that represent value type in first column.
        TabularDerivativeSmooth(const std::vector<FenestrationCommon::point> & values,
                                Variable property);

        //! \brief Construction of tabular derivative from standard vector values.
        //!
        //! \param list Initializer list used to construct tabular derivative.
        //! \param property Variable that represent value type in first column.
        TabularDerivativeSmooth(const std::initializer_list<FenestrationCommon::point> & list,
                                Variable property);

    protected:
        std::vector<FenestrationCommon::point> m_Curve;

        //! Overriden evaluation function.
        double evaluateFunction(double t_position, double t_previousTimestep) const override;

        //! Helper function that returns two closest points for interpolation.
        virtual std::pair<FenestrationCommon::point, FenestrationCommon::point>
          getInterpolationPoints(std::vector<FenestrationCommon::point>::const_iterator & it, const std::vector<FenestrationCommon::
                                 point> & table) const;
    };

    //////////////////////////////////////////////////////////////////
    ///  LiquidTransportationCurve
    //////////////////////////////////////////////////////////////////

    //! \brief Liquid transport coefficient curve K_l(w).
    //!
    //! Liquid transport curves are tabulated against material water content `w`.
    //! At a material interface, a node carries an averaged water content
    //! (computed across all materials touching that node by `Node2D::calcWaterContent`),
    //! which does NOT correspond to any single material's `w(phi)` at the node's
    //! humidity. Looking up `K_l` at the averaged `w` therefore yields a coefficient
    //! that is wrong for the element doing the lookup. To prevent that, the curve
    //! is bound to a specific material at construction time and resolves `w` for
    //! the lookup by calling `material.waterContent(node)` — i.e., the material's
    //! own `w(node.humidity)` — instead of going through the node's averaged
    //! `Variable::water` property. The bound material is normally the owning
    //! element's material.
    class LiquidTransportationCurve : public TabularFunction1D
    {
    public:
        //! Construction from standard vector values, bound to a material.
        //!
        //! \param vec      Liquid-transport curve points (w, K_l).
        //! \param material Material whose water content is used to resolve the
        //!                 lookup key when this curve is evaluated at a node.
        LiquidTransportationCurve(const std::vector<FenestrationCommon::point> & vec,
                                  const IMaterial & material);

        //! Construction from initializer list, bound to a material.
        //!
        //! \param list     Liquid-transport curve points (w, K_l).
        //! \param material Material whose water content is used to resolve the
        //!                 lookup key when this curve is evaluated at a node.
        LiquidTransportationCurve(const std::initializer_list<FenestrationCommon::point> & list,
                                  const IMaterial & material);

        //! Override that resolves `w` through the bound material instead of via
        //! `node.property(Variable::water)`. See class comment for the reason.
        double value(const INode2D & node) const override;

    protected:
        //! Helper function that returns two closest points for interpolation.
        std::pair<FenestrationCommon::point, FenestrationCommon::point> getInterpolationPoints(
            std::vector<FenestrationCommon::point>::const_iterator & it,
            const std::vector<FenestrationCommon::point> & table) const override;

    private:
        const IMaterial & m_Material;
    };

    //////////////////////////////////////////////////////////////////
    ///  SaturationFunction
    //////////////////////////////////////////////////////////////////

    //! \brief Simple saturation function.
    class SaturationFunction : public IFunction
    {
    public:
        //! Construction of saturation function.
        SaturationFunction();

    private:
        //! Overriden evaluation function.
        double evaluateFunction(double t_position, double t_previousTimestep) const override;
    };

    //////////////////////////////////////////////////////////////////
    ///  Heat of evaporation
    //////////////////////////////////////////////////////////////////

    //! \brief Heat of evaporation class represents heat evaporation dependency on temperature.
    class HeatOfEvaporation : public IFunction
    {
    public:
        //! Heat evaporation construction.
        HeatOfEvaporation();

    protected:
        //! Overriden evaluation function.
        double evaluateFunction(double t_position, double t_previousTimestep) const override;
    };


    //////////////////////////////////////////////////////////////////
    ///  Phase change
    //////////////////////////////////////////////////////////////////

    //! \brief Handles phase change around freezing point
    class PhaseChange : public IFunction
    {
    public:
        PhaseChange();

    protected:
        double evaluateFunction(double t_position, double t_previousTimestep) const override;
    };

}   // namespace HygroThermFEM
