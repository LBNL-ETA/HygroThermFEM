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
        virtual ~IValue() = default;

        virtual std::unique_ptr<IValue> clone() const = 0;
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
        IOperation(iValue t_Val1, iValue t_Val2, Operation t_Operation);

        double value(const State & state) const override;

    private:
        /// Do not want to make clone of operator publicly available
        virtual iValue clone() const override;

        /// Functions can be shared between different operations and that is why it is
        /// necessary to share function
        iValue m_Function1;
        iValue m_Function2;

        Operation m_Operation;

        /// This hold four basic operators (+, -. *. /) which is used to determine
        /// which function pointer is to be called
        std::map<Operation, std::function<double(double, double)>> m_Operator;
    };

    //////////////////////////////////////////////////////////////////
    ///  Operators
    //////////////////////////////////////////////////////////////////

    iValue operator+(iValue & left, iValue & right);

    iValue operator+(const double left, iValue & right);

    iValue operator+(iValue & left, const double right);

    iValue operator-(iValue & left, iValue & right);

    iValue operator-(const double left, iValue & right);

    iValue operator-(iValue & left, const double right);

    iValue operator*(iValue & left, iValue & right);

    iValue operator*(const double left, iValue & right);

    iValue operator*(iValue & left, const double right);

    iValue operator/(iValue & left, iValue & right);

    iValue operator/(const double left, iValue & right);

    iValue operator/(iValue & left, const double right);

    //////////////////////////////////////////////////////////////////
    ///  Constant
    //////////////////////////////////////////////////////////////////

    /// Simple constant curve.
    class Constant : public IFunction
    {
    public:
        static std::unique_ptr<Constant> create(const double value);

        virtual std::unique_ptr<IValue> clone() const override;

    private:
        Constant(const double value);

        double evaluateFunction(const double t_position) const override;

        double m_Value;
    };

	//////////////////////////////////////////////////////////////////
	///  State value
	//////////////////////////////////////////////////////////////////

	class StateValue : public IFunction
	{
	public:
		static std::unique_ptr<StateValue> create( Property property );

		virtual std::unique_ptr<IValue> clone() const override;
	private:
		StateValue(Property property);
		double evaluateFunction(const double t_position) const override;
	};

    //////////////////////////////////////////////////////////////////
    ///  Derivative
    //////////////////////////////////////////////////////////////////
    class Derivative : public IValue
    {
    public:
        static std::unique_ptr<Derivative> create(const iValue & t_Function);

        double value(const State & state) const override;

        virtual std::unique_ptr<IValue> clone() const override;

    private:
        Derivative(const iValue & t_Function);

        iValue m_Function;
    };

    //////////////////////////////////////////////////////////////////
    ///  TabularFunctions
    //////////////////////////////////////////////////////////////////

    /// Interface for classic tabular curve. There are different interpolation
    /// strategies and this is base class for all of them.
    class TabularFunction : public IFunction
    {
    public:
        static std::unique_ptr<TabularFunction>
          create(const std::vector<std::pair<double, double>> & values,
                 Property property,
                 const FenestrationCommon::Interpolator & interpolator =
                   FenestrationCommon::Interpolation::Linear);

        static std::unique_ptr<TabularFunction>
          create(const std::initializer_list<std::pair<double, double>> & list,
                 Property property,
                 const FenestrationCommon::Interpolator & interpolator =
                   FenestrationCommon::Interpolation::Linear);

        double max() const;

        double min() const;

        std::vector<std::pair<double, double>> getCurve() const;

        virtual std::unique_ptr<IValue> clone() const override;

    protected:
        TabularFunction(const std::vector<std::pair<double, double>> & values,
                        const Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Linear);

        TabularFunction(const std::initializer_list<std::pair<double, double>> & list,
                        const Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Linear);

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
        static std::unique_ptr<TabularDerivative>
          create(const std::vector<std::pair<double, double>> & values, Property property);
        static std::unique_ptr<TabularDerivative>
          create(const std::initializer_list<std::pair<double, double>> & list, Property property);

        virtual std::unique_ptr<IValue> clone() const override;

    protected:
        TabularDerivative(const std::vector<std::pair<double, double>> & values, Property property);

        TabularDerivative( const std::initializer_list< std::pair< double, double>> & list,
						   Property property );

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
        static std::unique_ptr<SuctionFunction>
          create(const std::vector<std::pair<double, double>> & values,
                 Property property,
                 const FenestrationCommon::Interpolator & interpolator =
                   FenestrationCommon::Interpolation::Logarithmic);

        static std::unique_ptr<SuctionFunction>
          create(const std::initializer_list<std::pair<double, double>> & list,
                 Property property,
                 const FenestrationCommon::Interpolator & interpolator =
                   FenestrationCommon::Interpolation::Logarithmic);

        virtual std::unique_ptr<IValue> clone() const override;

    protected:
        SuctionFunction(const std::vector<std::pair<double, double>> & values,
                        Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Logarithmic);

        SuctionFunction(const std::initializer_list<std::pair<double, double>> & list,
                        Property property,
                        const FenestrationCommon::Interpolator & interpolator =
                          FenestrationCommon::Interpolation::Logarithmic);

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
        static std::unique_ptr<SaturationFunction> create(Property property,
                                                          double saturationCoefficient = 9.2);

        virtual std::unique_ptr<IValue> clone() const override;

    private:
        SaturationFunction(Property property, double saturationCoefficient = 9.2);

        double evaluateFunction(const double t_position) const override;
        const double m_SaturationCoefficient;
    };

}   // namespace MoisThermFEM