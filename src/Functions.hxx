#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "Interpolator.hxx"

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
    };

    //////////////////////////////////////////////////////////////////
    ///  IFunction
    //////////////////////////////////////////////////////////////////

    /// Interface for functions
    class IFunction : public IValue
    {
    public:
        IFunction(Property t_Property);

        virtual double value(const State & state) const;

    protected:
        virtual double evaluateFunction(const double t_position = 0) const = 0;

        /// Property that is used to calculate function value. It is extracted from current
        /// domain (material) point.
        const Property m_Property;
    };

    //////////////////////////////////////////////////////////////////
    ///  IOperation
    //////////////////////////////////////////////////////////////////

    /// Entire class is used to mimic operator functions so that FEM functions can be
    /// written as ordinary equations.
    class IOperation : public IValue
    {
    public:
        IOperation(std::shared_ptr<IValue> & t_Val1,
                   std::shared_ptr<IValue> & t_Val2,
                   Operation t_Operation);

    public:
        double value(const State & state) const override;

    private:
        /// Functions can be shared between different operations and that is why it is
        /// necessary to share function
        std::shared_ptr< IValue > m_Function1;
        std::shared_ptr< IValue > m_Function2;

        Operation m_Operation;

        /// This hold four basic operators (+, -. *. /) which is used to determine
        /// which function pointer is to be called
        std::map<Operation, std::function<double(double, double)>> m_Operator;
    };

    //////////////////////////////////////////////////////////////////
    ///  Operators
    //////////////////////////////////////////////////////////////////

    std::shared_ptr<IValue> operator+(std::shared_ptr<IValue> & left,
                                      std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator+(const double left, std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator+(std::shared_ptr<IValue> & left, const double right);

    std::shared_ptr<IValue> operator-(std::shared_ptr<IValue> & left,
                                      std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator-(const double left, std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator-(std::shared_ptr<IValue> & left, const double right);

    std::shared_ptr<IValue> operator*(std::shared_ptr<IValue> & left,
                                      std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator*(const double left, std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator*(std::shared_ptr<IValue> & left, const double right);

    std::shared_ptr<IValue> operator/(std::shared_ptr<IValue> & left,
                                      std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator/(const double left, std::shared_ptr<IValue> & right);

    std::shared_ptr<IValue> operator/(std::shared_ptr<IValue> & left, const double right);

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
    ///  Derivative
    //////////////////////////////////////////////////////////////////
    class Derivative : public IValue
    {
    public:
        Derivative(std::shared_ptr<IValue> & t_Function);

        double value(const State & state) const override;

    private:
        std::shared_ptr<IValue> m_Function;
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
                        Property property,
                        FenestrationCommon::Interpolator interpolator =
                          FenestrationCommon::Interpolation::Linear);

        TabularFunction(std::initializer_list<std::pair<double, double>> & list,
                        Property property,
                        FenestrationCommon::Interpolator interpolator =
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

        TabularDerivative(std::initializer_list<std::pair<double, double>> & list,
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