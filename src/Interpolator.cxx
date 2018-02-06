#include <cmath>
#include <memory>

#include "Interpolator.hxx"

namespace FenestrationCommon {

	double IInterpolator::f( const std::pair< double, double > & t_point1,
													 const std::pair< double, double > & t_point2,
													 const double t_position ) const {
		return ( t_position - t_point1.first ) / ( t_point2.first - t_point1.first );
	}

	double LinearInt::interpolate( const std::pair< double, double > & t_point1,
																 const std::pair< double, double > & t_point2,
																 const double t_position ) const {
		const double f_pos = f( t_point1, t_point2, t_position );
		return f_pos * t_point2.second + ( 1 - f_pos ) * t_point1.second;
	}

	double LogarithmicInt::interpolate( const std::pair< double, double > & t_point1,
																			const std::pair< double, double > & t_point2,
																			const double t_position ) const {
		const double f_pos = f( t_point1, t_point2, t_position );
		return std::pow( t_point2.second, f_pos ) * std::pow( t_point1.second, ( 1 - f_pos ) );
	}

	std::unique_ptr< IInterpolator >
	InterpolatorFactory::getInterpolator( const Interpolator t_Interpolator ) {
		std::unique_ptr< IInterpolator > result = nullptr;
		switch ( t_Interpolator ) {
			case Interpolator::Linear:
				result = std::unique_ptr< IInterpolator >( new LinearInt() );
				break;
			case Interpolator::Logarithmic:
				result = std::unique_ptr< IInterpolator >( new LogarithmicInt() );
				break;
		}
		return result;
	}
}