#pragma once

#include <vector>
#include <map>
#include <functional>

#include "Interpolator.hxx"


/// Functions interface is used to build function that are used for matrix building. Functions are
/// stacked together to make full function that later will be stored in FEM element.

namespace MoisThermFEM {

	enum class Property;

	class State;

	enum class Operation {
		MULT, DIV, ADD, SUB
	};

	//////////////////////////////////////////////////////////////////
	///  IFunction
	//////////////////////////////////////////////////////////////////

	/// Interface for functions
	class IFunction {
	public:
		IFunction( Property t_Property );

		virtual double value( const State & state ) const = 0;

	protected:
		virtual double getValue( const double t_position ) const = 0;

		Property m_Property;

	};

	//////////////////////////////////////////////////////////////////
	///  IDecoratingFunction
	//////////////////////////////////////////////////////////////////

	class IDecoratingFunction : public IFunction {
	public:

		IDecoratingFunction( Property property );

		IDecoratingFunction( Property property, std::unique_ptr< IFunction > & m_Curve,
												 Operation operation = Operation::MULT );

		double value( const State & state ) const final;

	protected:
		std::unique_ptr< IFunction > m_Function;
		Operation m_Operation;
		std::map< Operation, std::function< double( double, double ) > > m_Operator;

	};

	//////////////////////////////////////////////////////////////////
	///  Constant
	//////////////////////////////////////////////////////////////////

	/// Simple constant curve.
	class Constant : public IDecoratingFunction {
	public:
		Constant( const double value, Property property );

		Constant( const double value, Property property, std::unique_ptr< IFunction > & t_Curve,
							Operation operation = Operation::MULT );

	private:
		double getValue( const double t_position = 0 ) const override;

		double m_Value;
	};

	//////////////////////////////////////////////////////////////////
	///  TabularFunctions
	//////////////////////////////////////////////////////////////////

	/// Interface for classic tabular curve. There are different interpolation strategies and
	/// this is base class for all of them.
	class TabularFunction : public IDecoratingFunction {
	public:
		TabularFunction( const std::vector< std::pair< double, double > > & values,
										 Property property,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		TabularFunction( std::initializer_list< std::pair< double, double>> & list,
										 Property property,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		TabularFunction( const std::vector< std::pair< double, double > > & values,
										 Property property,
										 std::unique_ptr< IFunction > & t_Curve,
										 Operation operation = Operation::MULT,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		TabularFunction( const std::initializer_list< std::pair< double, double > > & list,
										 Property property,
										 std::unique_ptr< IFunction > & t_Curve,
										 Operation operation = Operation::MULT,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		double max() const;

		double min() const;

	protected:
		std::vector< std::pair< double, double > > m_Curve;
		FenestrationCommon::Interpolator m_Interpolator;

		double getValue( const double t_position ) const override;

		virtual std::pair< std::pair< double, double >, std::pair< double, double > >
		getInterpolationPoints(
				std::vector< std::pair< double, double > >::const_iterator & it ) const;
	};

	//////////////////////////////////////////////////////////////////
	///  SuctionFunction
	//////////////////////////////////////////////////////////////////

	/// Class that behaves like suction curve. It is standard (linear or logarithmic) interpolation
	/// except for the results in first range where curve will return constant value equal to the
	/// first point
	class SuctionFunction : public TabularFunction {
	public:
		SuctionFunction( const std::vector< std::pair< double, double > > & values,
										 Property property,
										 const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

		SuctionFunction( const std::initializer_list< std::pair< double, double > > & list,
										 Property property,
										 const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

		SuctionFunction( const std::vector< std::pair< double, double > > & values,
										 Property property,
										 std::unique_ptr< IFunction > & t_Curve,
										 Operation operation = Operation::MULT,
										 const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

		SuctionFunction( const std::initializer_list< std::pair< double, double > > & list,
										 Property property,
										 std::unique_ptr< IFunction > & t_Curve,
										 Operation operation = Operation::MULT,
										 const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

	protected:
		std::pair< std::pair< double, double >, std::pair< double, double > >
		getInterpolationPoints(
				std::vector< std::pair< double, double > >::const_iterator & it ) const override;

	};

	//////////////////////////////////////////////////////////////////
	///  FirstDerivativeFunction
	//////////////////////////////////////////////////////////////////
	class FirstDerivativeFunction : public TabularFunction {

	public:
		FirstDerivativeFunction( const std::vector< std::pair< double, double > > & values,
														 Property property,
														 const FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		FirstDerivativeFunction( const std::initializer_list< std::pair< double, double > > & list,
														 Property property,
														 const FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		FirstDerivativeFunction( const std::vector< std::pair< double, double > > & values,
														 Property property,
														 std::unique_ptr< IFunction > & t_Curve,
														 Operation operation = Operation::MULT,
														 const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Linear );

		FirstDerivativeFunction( const std::initializer_list< std::pair< double, double > > & list,
														 Property property,
														 std::unique_ptr< IFunction > & t_Curve,
														 Operation operation = Operation::MULT,
														 const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Linear );

	private:
		double getValue( const double t_position ) const override;

		virtual double firstDerivative( const double t_position ) const final;

	};

	//////////////////////////////////////////////////////////////////
	///  SaturationFunction
	//////////////////////////////////////////////////////////////////

	/// Simple constant curve.
	class SaturationFunction : public IDecoratingFunction {
	public:
		SaturationFunction( Property property );

		SaturationFunction( Property property, std::unique_ptr< IFunction > & t_Curve,
														Operation operation );

	private:
		double getValue( const double t_position ) const override;

	};

}