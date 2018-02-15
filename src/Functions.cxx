#include <algorithm>
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
																						std::unique_ptr< IFunction > & m_Curve )
			: IFunction( property ), m_Function( std::move( m_Curve ) ) {}

	double IDecoratingFunction::value( const State & state ) const {
		auto value = state.getValue( m_Property );
		return m_Function != nullptr ? getValue( value ) * m_Function->value( state )
																 : getValue( value );
	}

	//////////////////////////////////////////////////////////////////
	///  Constant
	//////////////////////////////////////////////////////////////////

	Constant::Constant( const double value, Property property ) : IDecoratingFunction( property ), m_Value( value ) {}

	double Constant::getValue( const double ) const {
		return m_Value;
	}

	Constant::Constant( const double value, Property property, std::unique_ptr< IFunction > & t_Curve )
			: IDecoratingFunction( property, t_Curve ), m_Value( value ) {

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
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			FenestrationCommon::Interpolator interpolator ) : IDecoratingFunction( property ), m_Curve( list ),
																												m_Interpolator(
																														std::move( interpolator ) ) {}

	TabularFunction::TabularFunction( const std::vector< std::pair< double, double > > & values,
																		Property property,
																		std::unique_ptr< IFunction > & t_Curve,
																		FenestrationCommon::Interpolator interpolator )
			: IDecoratingFunction( property, t_Curve ), m_Curve( values ),
				m_Interpolator( std::move( interpolator ) ) {

	}

	TabularFunction::TabularFunction(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			std::unique_ptr< IFunction > & t_Curve, FenestrationCommon::Interpolator interpolator ) :
			IDecoratingFunction( property, t_Curve ), m_Curve( list ), m_Interpolator( std::move( interpolator ) ) {

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

	//////////////////////////////////////////////////////////////////
	///  SuctionCurve
	//////////////////////////////////////////////////////////////////

	SuctionCurve::SuctionCurve( const std::vector< std::pair< double, double > > & values,
															Property property,
															const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( values, property, interpolator ) {

	}

	SuctionCurve::SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
															Property property,
															const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( list, property, interpolator ) {

	}

	SuctionCurve::SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
															Property property,
															std::unique_ptr< IFunction > & t_Curve,
															const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( list, property, t_Curve, interpolator ) {

	}

	SuctionCurve::SuctionCurve( const std::vector< std::pair< double, double > > & values,
															Property property,
															std::unique_ptr< IFunction > & t_Curve,
															const FenestrationCommon::Interpolator & interpolator )
			: TabularFunction( values, property, t_Curve, interpolator ) {

	}

	std::pair< std::pair< double, double >, std::pair< double, double>>
	SuctionCurve::getInterpolationPoints(
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

	FirstDerivativeCurve::FirstDerivativeCurve(
			const std::vector< std::pair< double, double > > & values,
			Property property,
			const FenestrationCommon::Interpolator interpolator )
			: TabularFunction( values, property, interpolator ) {

	}

	FirstDerivativeCurve::FirstDerivativeCurve(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			const FenestrationCommon::Interpolator interpolator ) :
			TabularFunction( list, property, interpolator ) {

	}

	FirstDerivativeCurve::FirstDerivativeCurve(
			const std::vector< std::pair< double, double > > & values,
			Property property,
			std::unique_ptr< IFunction > & t_Curve, const FenestrationCommon::Interpolator interpolator )
			:
			TabularFunction( values, property, t_Curve, interpolator ) {

	}

	FirstDerivativeCurve::FirstDerivativeCurve(
			const std::initializer_list< std::pair< double, double > > & list,
			Property property,
			std::unique_ptr< IFunction > & t_Curve, const FenestrationCommon::Interpolator interpolator )
			:
			TabularFunction( list, property, t_Curve, interpolator ) {

	}

	double FirstDerivativeCurve::getValue( const double t_position ) const {
		return firstDerivative( t_position );
	}

	double FirstDerivativeCurve::firstDerivative( const double t_position ) const {
		const double small = 1e-8;
		double val1 = TabularFunction::getValue( t_position );
		double val2 = TabularFunction::getValue( t_position + small );
		return ( val2 - val1 ) / small;
	}

}

