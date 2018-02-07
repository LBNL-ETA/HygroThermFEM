#pragma once

#include <utility>
#include <map>
#include <functional>

namespace FenestrationCommon {

	enum class Interpolation {
		Linear, Logarithmic
	};

	class Interpolator {
	public:
		Interpolator( const Interpolation t_interpolation );

		double interpolate( const std::pair< double, double > & t_point1,
												const std::pair< double, double > & t_point2,
												const double t_position ) const;

	private:
		double f( const std::pair< double, double > & t_point1,
							const std::pair< double, double > & t_point2,
							const double t_position ) const;

		std::map< Interpolation, std::function< double( const std::pair< double, double > & t_point1,
																									 const std::pair< double, double > & t_point2,
																									 const double t_position ) > > m_Functions;

		Interpolation m_Interpolation;
	};


}