#pragma once

#include <vector>

#include "Interpolator.hxx"

namespace FenestrationCommon {

	//////////////////////////////////////////////////////////////////
	///  IFunction
	//////////////////////////////////////////////////////////////////

	/// Interface for curves.
	class IFunction {
	public:
		virtual double value( const double t_position = 0 ) const = 0;

	protected:
		virtual double getValue( const double t_position ) const = 0;

	};

	//////////////////////////////////////////////////////////////////
	///  IDecoratingFunction
	//////////////////////////////////////////////////////////////////

	class IDecoratingFunction : public IFunction {
	public:

		IDecoratingFunction();

		IDecoratingFunction( std::unique_ptr< IFunction > & m_Curve );

		double value( const double t_position ) const final;

	protected:
		std::unique_ptr< IFunction > m_Function;
	};

	//////////////////////////////////////////////////////////////////
	///  Constant
	//////////////////////////////////////////////////////////////////

	/// Simple constant curve.
	class Constant : public IDecoratingFunction {
	public:
		Constant( const double value );
		Constant( const double value, std::unique_ptr< IFunction > & t_Curve );

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
					 Interpolator interpolator = Interpolation::Linear );

		TabularFunction( const std::initializer_list< std::pair< double, double > > & list,
					 Interpolator interpolator = Interpolation::Linear );

		TabularFunction( const std::vector< std::pair< double, double > > & values,
					 std::unique_ptr< IFunction > & t_Curve,
					 Interpolator interpolator = Interpolation::Linear );

		TabularFunction( const std::initializer_list< std::pair< double, double > > & list,
					 std::unique_ptr< IFunction > & t_Curve,
					 Interpolator interpolator = Interpolation::Linear );

	protected:
		std::vector< std::pair< double, double > > m_Curve;
		Interpolator m_Interpolator;

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
									const Interpolator & interpolator = Interpolation::Logarithmic );

		SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
									const Interpolator & interpolator = Interpolation::Logarithmic );

		SuctionCurve( const std::vector< std::pair< double, double > > & values,
									std::unique_ptr< IFunction > & t_Curve,
									const Interpolator & interpolator = Interpolation::Logarithmic );

		SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
									std::unique_ptr< IFunction > & t_Curve,
									const Interpolator & interpolator = Interpolation::Logarithmic );

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
					 const Interpolator interpolator = Interpolation::Linear );

		FirstDerivativeCurve( const std::initializer_list< std::pair< double, double > > & list,
					 const Interpolator interpolator = Interpolation::Linear );

		FirstDerivativeCurve( const std::vector< std::pair< double, double > > & values,
					 std::unique_ptr< IFunction > & t_Curve,
					 const Interpolator interpolator = Interpolation::Linear );

		FirstDerivativeCurve( const std::initializer_list< std::pair< double, double > > & list,
					 std::unique_ptr< IFunction > & t_Curve,
					 const Interpolator interpolator = Interpolation::Linear );

	private:
		double getValue( const double t_position ) const override;
		virtual double firstDerivative( const double t_position ) const final;

	};
}