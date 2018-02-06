#pragma once

#include <utility>

namespace FenestrationCommon {

	enum class Interpolator { Linear, Logarithmic };

	class IInterpolator {
	public:
		IInterpolator() = default;

		virtual double interpolate( const std::pair< double, double > & t_point1,
						const std::pair< double, double > & t_point2,
						const double t_position ) const = 0;

	protected:
		double f(const std::pair< double, double > & t_point1,
						 const std::pair< double, double > & t_point2,
						 const double t_position) const;
	};

	class LinearInt : public IInterpolator {
	public:
		LinearInt() = default;

		double interpolate( const std::pair< double, double > & t_point1,
												const std::pair< double, double > & t_point2,
												const double t_position ) const override;

	};

	class LogarithmicInt : public IInterpolator {
	public:
		LogarithmicInt() = default;

		double interpolate( const std::pair< double, double > & t_point1,
												const std::pair< double, double > & t_point2,
												const double t_position ) const override;
	};

	class InterpolatorFactory {
	public:

		static std::unique_ptr< IInterpolator > getInterpolator( const Interpolator t_Interpolator );


	};


}