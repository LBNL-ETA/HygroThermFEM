#include <algorithm>
#include <cmath>

#include "Functions.hxx"
#include "State.hxx"

namespace MoisThermFEM {

	//////////////////////////////////////////////////////////////////
	///  IDecoratingCurve
	//////////////////////////////////////////////////////////////////

	IFunction::IFunction( Property t_Property ) : m_Property( t_Property ) {}

	//////////////////////////////////////////////////////////////////
	///  IDecoratingCurve
	//////////////////////////////////////////////////////////////////

	IDecoratingFunction::IDecoratingFunction( Property property ) : IFunction( property ),
																																	m_Function{ nullptr } {}

	IDecoratingFunction::IDecoratingFunction( Property property,
																						std::unique_ptr< IFunction > & m_Curve,
																						Operation operation )
			: IFunction( property ), m_Function( std::move( m_Curve ) ), m_Operation( operation ) {
		m_Operator[ Operation::MULT ] = [ & ]( double a, double b ) { return a * b; };
		m_Operator[ Operation::DIV ] = [ & ]( double a, double b ) { return a / b; };
		m_Operator[ Operation::ADD ] = [ & ]( double a, double b ) { return a + b; };
		m_Operator[ Operation::SUB ] = [ & ]( double a, double b ) { return a - b; };
	}

	double IDecoratingFunction::value( const State & state ) const {
		auto value = state.getValue( m_Property );
		return m_Function != nullptr ? m_Operator.at( m_Operation )( getValue( value ),
																																 m_Function->value( state ) )
																 : getValue( value );
	}

	//////////////////////////////////////////////////////////////////
	///  Constant
	//////////////////////////////////////////////////////////////////

	Constant::Constant( const double value, Property property ) : IDecoratingFunction( property ),
																																m_Value( value ) {}

	double Constant::getValue( const double ) const {
		return m_Value;
	}

	Constant::Constant( const double value, Property property,
											std::unique_ptr< IFunction > & t_Curve,
											Operation operation )
			: IDecoratingFunction( property, t_Curve, operation ), m_Value( value ) {

	}

	//////////////////////////////////////////////////////////////////
	///  Curve
	//////////////////////////////////////////////////////////////////

	TabularFunction::TabularFunction( const std::vector< std::pair< double, double > > & values,
																		Property property,
																		FenestrationCommon::Interpolator interpolator )
			: IDecoratingFunction( property ), m_Curve( values ),
				m_Interpolator( std::move( interpolator ) ) {}

	TabularFunction::TabularFunction(
			std::initializer_list< std::pair< double, double > > & list,
			Property property,
			FenestrationCommon::Interpolator interpolator ) : IDecoratingFunction( property ),
																												m_Curve( std::move( list ) ),
																												m_Interpolator(
																														std::move( interpolator ) ) {}

	TabularFunction::TabularFunction( const std::vector< std::pair< double, double > > & values,
																		Property property,
																		std::unique_ptr< IFunction > & t_Curve,
																		Operation operation,
																		FenestrationCommon::Interpolator interpolator )
			: IDecoratingFunction( property, t_Curve, operation ), m_Curve( values ),
				m_Interpolator( std::move( interpolator ) ) {

	}

	TabularFunction::TabularFunction(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			std::unique_ptr< IFunction > & t_Curve,
			Operation operation,
			FenestrationCommon::Interpolator interpolator ) :
			IDecoratingFunction( property, t_Curve, operation ), m_Curve( list ),
			m_Interpolator( std::move( interpolator ) ) {

	}

	double TabularFunction::getValue( const double t_position ) const {
		auto it = std::find_if( m_Curve.begin(), m_Curve.end(),
														[ & ]( std::pair< double, double > val ) {
															return val.first > t_position;
														} );

		auto points = getInterpolationPoints( it );

		return m_Interpolator.interpolate( points.first, points.second, t_position );
	}

	std::pair< std::pair< double, double >, std::pair< double, double >>
	TabularFunction::getInterpolationPoints(
			std::vector< std::pair< double, double > >::const_iterator & it ) const {
		auto pt2 = it == m_Curve.end() ? m_Curve.back() : *it;
		if( it != m_Curve.begin() ) {
			--it;
		}
		auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

		return std::make_pair( pt1, pt2 );
	}

	double TabularFunction::max() const {
		return m_Curve.back().second;
	}

	double TabularFunction::min() const {
		return m_Curve.front().second;
	}

	//////////////////////////////////////////////////////////////////
	///  SuctionCurve
	//////////////////////////////////////////////////////////////////

	SuctionFunction::SuctionFunction( const std::vector< std::pair< double, double > > & values,
																		Property property,
																		const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( values, property, interpolator ) {

	}

	SuctionFunction::SuctionFunction(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( list, property, interpolator ) {

	}

	SuctionFunction::SuctionFunction(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			std::unique_ptr< IFunction > & t_Curve,
			Operation operation,
			const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( list, property, t_Curve, operation, interpolator ) {

	}

	SuctionFunction::SuctionFunction( const std::vector< std::pair< double, double > > & values,
																		Property property,
																		std::unique_ptr< IFunction > & t_Curve,
																		Operation operation,
																		const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( values, property, t_Curve, operation, interpolator ) {

	}

	std::pair< std::pair< double, double >, std::pair< double, double > >
	SuctionFunction::getInterpolationPoints(
			std::vector< std::pair< double, double > >::const_iterator & it ) const {

		/// Suction curve takes care that first segment of curve always return value of first element.
		it == m_Curve.end() ? m_Curve.back() : *it;
		auto second = m_Curve.begin() + 1;
		auto pt2 = it == second ? m_Curve.front() : *it;
		if( it != m_Curve.begin() ) {
			--it;
		}

		auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

		return std::make_pair( pt1, pt2 );
	}

	//////////////////////////////////////////////////////////////////
	///  FirstDerivativeCurve
	//////////////////////////////////////////////////////////////////

	FirstDerivativeFunction::FirstDerivativeFunction(
			const std::vector< std::pair< double, double > > & values,
			Property property,
			const FenestrationCommon::Interpolator interpolator )
			: TabularFunction( values, property, interpolator ) {

	}

	FirstDerivativeFunction::FirstDerivativeFunction(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			const FenestrationCommon::Interpolator interpolator ) :
			TabularFunction( list, property, interpolator ) {

	}

	FirstDerivativeFunction::FirstDerivativeFunction(
			const std::vector< std::pair< double, double > > & values,
			Property property,
			std::unique_ptr< IFunction > & t_Curve, Operation operation,
			const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( values, property, t_Curve, operation, interpolator ) {

	}

	FirstDerivativeFunction::FirstDerivativeFunction(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			std::unique_ptr< IFunction > & t_Curve,
			Operation operation,
			const FenestrationCommon::Interpolator & interpolator ) :
			TabularFunction( list, property, t_Curve, operation, interpolator ) {

	}

	double FirstDerivativeFunction::getValue( const double t_position ) const {
		return firstDerivative( t_position );
	}

	double FirstDerivativeFunction::firstDerivative( const double t_position ) const {
		const double small = 1e-8;
		double val1 = TabularFunction::getValue( t_position );
		double val2 = TabularFunction::getValue( t_position + small );
		return ( val2 - val1 ) / small;
	}

	//////////////////////////////////////////////////////////////////
	///  SaturationFunction
	//////////////////////////////////////////////////////////////////

	SaturationFunction::SaturationFunction( Property property ) :
			IDecoratingFunction( property ) {}

	SaturationFunction::SaturationFunction( Property property, std::unique_ptr< IFunction > & t_Curve,
																					Operation operation )
			: IDecoratingFunction( property, t_Curve, operation ) {

	}

	double SaturationFunction::getValue( const double t_position ) const {
		auto temp = 77.345 + 0.0057 * t_position - 7235 / t_position;
		temp = std::exp( temp );
		return temp / ( 461.4 * std::pow( t_position, 9.2 ) );
	}

}

