#pragma once

#include <vector>

#include "Interpolator.hxx"

namespace MoisThermFEM {

	enum class Property;
	class State;

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

		IDecoratingFunction( Property property, std::unique_ptr< IFunction > & m_Curve );

		double value( const State & state ) const final;

	protected:
		std::unique_ptr< IFunction > m_Function;
	};

	//////////////////////////////////////////////////////////////////
	///  Constant
	//////////////////////////////////////////////////////////////////

	/// Simple constant curve.
	class Constant : public IDecoratingFunction {
	public:
		Constant( const double value, Property property );

		Constant( const double value, Property property, std::unique_ptr< IFunction > & t_Curve );

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

		TabularFunction( const std::initializer_list< std::pair< double, double > > & list,
										 Property property,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		TabularFunction( const std::vector< std::pair< double, double > > & values,
										 Property property,
										 std::unique_ptr< IFunction > & t_Curve,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		TabularFunction( const std::initializer_list< std::pair< double, double > > & list,
										 Property property,
										 std::unique_ptr< IFunction > & t_Curve,
										 FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

	protected:
		std::vector< std::pair< double, double > > m_Curve;
		FenestrationCommon::Interpolator m_Interpolator;

		double getValue( const double t_position ) const override;

		virtual std::pair< std::pair< double, double >, std::pair< double, double > >
		getInterpolationPoints(
				std::vector< std::pair< double, double > >::const_iterator & it ) const;
	};

	//////////////////////////////////////////////////////////////////
	///  SuctionCurve
	//////////////////////////////////////////////////////////////////

	/// Class that behaves like suction curve. It is standard (linear or logarithmic) interpolation
	/// except for the results in first range where curve will return constant value equal to the
	/// first point
	class SuctionCurve : public TabularFunction {
	public:
		SuctionCurve( const std::vector< std::pair< double, double > > & values,
									Property property,
									const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

		SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
									Property property,
									const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

		SuctionCurve( const std::vector< std::pair< double, double > > & values,
									Property property,
									std::unique_ptr< IFunction > & t_Curve,
									const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

		SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
									Property property,
									std::unique_ptr< IFunction > & t_Curve,
									const FenestrationCommon::Interpolator & interpolator = FenestrationCommon::Interpolation::Logarithmic );

	protected:
		std::pair< std::pair< double, double >, std::pair< double, double > >
		getInterpolationPoints(
				std::vector< std::pair< double, double > >::const_iterator & it ) const override;

	};

	//////////////////////////////////////////////////////////////////
	///  FirstDerivativeCurve
	//////////////////////////////////////////////////////////////////
	class FirstDerivativeCurve : public TabularFunction {

	public:
		FirstDerivativeCurve( const std::vector< std::pair< double, double > > & values,
													Property property,
													const FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		FirstDerivativeCurve( const std::initializer_list< std::pair< double, double > > & list,
													Property property,
													const FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		FirstDerivativeCurve( const std::vector< std::pair< double, double > > & values,
													Property property,
													std::unique_ptr< IFunction > & t_Curve,
													const FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

		FirstDerivativeCurve( const std::initializer_list< std::pair< double, double > > & list,
													Property property,
													std::unique_ptr< IFunction > & t_Curve,
													const FenestrationCommon::Interpolator interpolator = FenestrationCommon::Interpolation::Linear );

	private:
		double getValue( const double t_position ) const override;

		virtual double firstDerivative( const double t_position ) const final;

	};
}