#include "Curve.hxx"

namespace FenestrationCommon {

	double ICurve::firstDerivative( const double t_position ) const {
		const double small = 1e-8;
		double val1 = value( t_position );
		double val2 = value( t_position + small );
		return ( val2 - val1 ) / small;
	};

	Curve::Curve( const std::vector< std::pair< double, double > > & values,
								const Interpolator interpolator ) : m_Curve( values ),
																										m_Interpolator( interpolator ) {}

	Curve::Curve( const std::initializer_list< std::pair< double, double > > & list,
								const Interpolator interpolator ) : m_Curve( list ),
																										m_Interpolator( interpolator ) {}

	double Curve::value( const double t_position ) const {
		auto it = std::find_if( m_Curve.begin(), m_Curve.end(),
														[ & ]( std::pair< double, double > val ) {
															return
																	val.first >
																	t_position;
														} );

		auto points = getInterpolationPoints( it );

		return m_Interpolator.interpolate( points.first, points.second, t_position );
	}

	std::pair< std::pair< double, double >, std::pair< double, double >>
	Curve::getInterpolationPoints(
			std::vector< std::pair< double, double > >::const_iterator & it ) const {
		auto pt2 = it == m_Curve.end() ? m_Curve.back() : *it;
		if( it != m_Curve.begin() ) {
			--it;
		}
		auto pt1 = it == m_Curve.begin() ? m_Curve.front() : *it;

		return std::make_pair( pt1, pt2 );
	}

	Constant::Constant( const double value ) : ICurve(), m_Value( value ) {}

	double Constant::value( const double ) const {
		return m_Value;
	}

	SuctionCurve::SuctionCurve( const std::vector< std::pair< double, double>> & values,
															const Interpolator & interpolator ) : Curve( values, interpolator ) {}

	SuctionCurve::SuctionCurve( const std::initializer_list< std::pair< double, double>> & list,
															const Interpolator & interpolator ) : Curve( list, interpolator ) {}

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
}

