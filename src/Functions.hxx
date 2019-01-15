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
	enum class Variable;

	class State;

	class Node2D;

	class INodes;

	//! Saturation function that will be used at boundary conditions
	double boundarySaturationAtTemperature(double temperature, double exponent = 8.2);

	//! Saturation function
	double saturationAtTemperature(double temperature, double exponent = 9.2);

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
	//! Passed Variable will be used to automatically read value of state variable and perform
	//! calculations.
	class IFunction : public IValue
	{
	public:
		//! Basic constructor
		IFunction(Variable t_Property   //!< Variable for which function will be calculated.
		);

		//! Returns function evaluation for given node.
		double value(const Node2D & node   //!< Node at which function will be evaluated.
		) const override;

	protected:
		//! Interface for function definition. This is place where in inherited classes function
		//! definitions will be stored.
		virtual double
		evaluateFunction(double t_position = 0,   //!< Value at which function will be evaluated.
		double t_previousTimestep = 0 //!< Value from previous timestep
		) const = 0;

		/// Variable that is used to calculate function value. It is extracted from current
		/// domain (material) point.
		const Variable m_Property;
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
		double evaluateFunction(double t_position, double t_previousTimestep) const override;

		double m_Value;
	};

	enum class Operation
	{
		MULT, DIV, ADD, SUB
	};

	//////////////////////////////////////////////////////////////////
	///  IOperation
	//////////////////////////////////////////////////////////////////

	//! \brief It is used to support operation with functions.
	//!
	//! Functions that are child of IValue are used directly in differential equations.
	//! Those functions can be stacked with ordinary operations. This class is used to support
	//! those operations.
	template <class T, class U>
	class IOperation : public IValue
	{
	public:
		//! Constructor that accept two operands and operation.
		IOperation(const T t, const U s, const Operation & op) : m_Function1(std::move(t)),
																 m_Function2(std::move(s)),
																 m_Operation(op)
		{
			m_Operator[Operation::MULT] = [&](double a, double b)
			{ return a * b; };
			m_Operator[Operation::DIV] = [&](double a, double b)
			{ return a / b; };
			m_Operator[Operation::ADD] = [&](double a, double b)
			{ return a + b; };
			m_Operator[Operation::SUB] = [&](double a, double b)
			{ return a - b; };
		}

		//! Returns value of operation.
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

	//! template operator + for any class derived from IValue. It accepts two operands.
	template <typename T, typename U>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
	operator+(const T & t, const U & u)
	{
		return IOperation<T, U>(t, u, Operation::ADD);
	}

	//! template operator + for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename T>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
	operator+(const T & t, const double & u)
	{
		Constant con{u};
		return IOperation<T, Constant>(t, con, Operation::ADD);
	}

	//! template operator + for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename U>
	typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
	operator+(const double & t, const U & u)
	{
		Constant con{t};
		return IOperation<Constant, U>(con, u, Operation::ADD);
	}

	////	operator-

	//! template operator - for any class derived from IValue. It accepts two operands.
	template <typename T, typename U>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
	operator-(const T & t, const U & u)
	{
		return IOperation<T, U>(t, u, Operation::SUB);
	}

	//! template operator - for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename T>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
	operator-(const T & t, const double & u)
	{
		Constant con{u};
		return IOperation<T, Constant>(t, con, Operation::SUB);
	}

	//! template operator - for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename U>
	typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
	operator-(const double & t, const U & u)
	{
		Constant con{t};
		return IOperation<Constant, U>(con, u, Operation::SUB);
	}

	////	operator*

	//! template operator * for any class derived from IValue. It accepts two operands.
	template <typename T, typename U>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
	operator*(const T & t, const U & u)
	{
		return IOperation<T, U>(t, u, Operation::MULT);
	}

	//! template operator * for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename T>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
	operator*(const T & t, const double & u)
	{
		Constant con{u};
		return IOperation<T, Constant>(t, con, Operation::MULT);
	}

	//! template operator * for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename U>
	typename std::enable_if<std::is_base_of<IValue, U>::value, IOperation<Constant, U>>::type
	operator*(const double & t, const U & u)
	{
		Constant con{t};
		return IOperation<Constant, U>(con, u, Operation::MULT);
	}

	////	operator/

	//! template operator / for any class derived from IValue. It accepts two operands.
	template <typename T, typename U>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, U>>::type
	operator/(const T & t, const U & u)
	{
		return IOperation<T, U>(t, u, Operation::DIV);
	}

	//! template operator / for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename T>
	typename std::enable_if<std::is_base_of<IValue, T>::value, IOperation<T, Constant>>::type
	operator/(const T & t, const double & u)
	{
		Constant con{u};
		return IOperation<T, Constant>(t, con, Operation::DIV);
	}

	//! template operator / for any class derived from IValue. It accepts two operands one of which
	//! is type of double.
	template <typename U>
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
		//! Constructor
		StateValue(Variable property   //!< Variable that StateValue represents.
		);

	private:
		//! Inherited function evaluation for current property
		double evaluateFunction(double t_position, double t_previousTimestep) const override;
	};

	//////////////////////////////////////////////////////////////////
	///  TabularFunction
	//////////////////////////////////////////////////////////////////

	//! \brief Interface for classic tabular curve.
	//!
	//! Tabular function is made to simply return adequate value from table. Search function is
	//! performed over first column and return value is calculated from second column.
	class TabularFunction : public IFunction
	{
	public:
		//! Table construction from standard vector.
		TabularFunction(
			const std::vector<std::pair<double, double>> & values,            //!< Vector of x,y pairs that will go into table.
			Variable property,   //!< Variable that represent value type in first column.
			const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Linear   //!< Interpolation strategy used for in
			//!< between values.
		);

		//! Table construction from initializer list.
		TabularFunction(
			const std::initializer_list<std::pair<double, double>> & list,              //!< Initializer list of x,y pairs that will go into table.
			Variable property,   //!< Variable that represent value type in first column.
			const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Linear   //!< Interpolation strategy used for in
			//!< between values.
		);

		//! Returns maximum value for first column.
		double maxX() const;

		//! Returns maximum value from second column.
		double maxY() const;

		//! Returns minimum value for first column.
		double minX() const;

		//! Returns minimum value from second column.
		double minY() const;

		//! Returns tabular values as standard vector of pairs.
		const std::vector<std::pair<double, double>> & getCurve() const;

	protected:
		std::vector<std::pair<double, double>> m_Curve;
		FenestrationCommon::Interpolator m_Interpolator;

		//! Overriden evaluation function.
		double evaluateFunction(double t_position, double t_previousTimestep) const override;

		//! Helper function that returns two closest points for interpolation.
		virtual std::pair<std::pair<double, double>, std::pair<double, double>>
		getInterpolationPoints(std::vector<std::pair<double, double>>::const_iterator & it) const;
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
		//! Construction of tabular derivative from standard vector values.
		TabularDerivative(
			const std::vector<std::pair<double, double>> & values, //!< Pair of vector values used to construct tabular derivative.
			Variable property //!< Variable that represent value type in first column.
		);

		//! Construction of tabular derivative from standard vector values.
		TabularDerivative(
			const std::initializer_list<std::pair<double, double>> & list, //!< Initializer list used to construct tabular derivative.
			Variable property //!< Variable that represent value type in first column.
		);

	protected:
		std::vector<std::pair<double, double>> m_Curve;

		//! Overriden evaluation function.
		double evaluateFunction(double t_position, double t_previousTimestep) const override;

		//! Helper function that returns two closest points for interpolation.
		virtual std::pair<std::pair<double, double>, std::pair<double, double>>
		getInterpolationPoints(std::vector<std::pair<double, double>>::const_iterator & it) const;
	};

	//////////////////////////////////////////////////////////////////
	///  SuctionFunction
	//////////////////////////////////////////////////////////////////

	//! \brief Sorption curve is specialized type for tabular function.
	//!
	//! Sorption curve is table that represents material water content as function of
	//! relative humidity. For in between values, logarithmic interpolation is used.
	class SuctionCurve : public TabularFunction
	{
	public:
		//! Construction of suction curve from standard vector values.
		SuctionCurve(
			const std::vector<std::pair<double, double>> & values //!< Sorption curve values in standard vector form.
		);

		//! Construction of suction curve from initializer list.
		SuctionCurve(
			const std::initializer_list<std::pair<double, double>> & list //!< Soprtion curve values in initializer list form.
		);

	protected:
		//! Helper function that returns two closest points for interpolation.
		std::pair<std::pair<double, double>, std::pair<double, double>> getInterpolationPoints(
			std::vector<std::pair<double, double>>::const_iterator & it) const override;
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

}   // namespace MoisThermFEM
