#pragma once

#include <vector>

#include "Interpolator.hxx"

namespace FenestrationCommon {

	/// Interface for curves. Satisfies basic need for value and first derivative.
	class ICurve {
	public:
		virtual double value( const double t_position ) const = 0;

		virtual double firstDerivative( const double t_position ) const final;
	};

	/// Interface for classic tabular curve. There are different interpolation strategies and
	/// this is base class for all of them.
	class Curve : public ICurve {
	public:
		Curve( const std::vector< std::pair< double, double > > & values,
					 const Interpolator interpolator = Interpolation::Linear );

		Curve( const std::initializer_list< std::pair< double, double > > & list,
					 const Interpolator interpolator = Interpolation::Linear );

		double value( const double t_position ) const override;

	protected:
		std::vector< std::pair< double, double > > m_Curve;
		Interpolator m_Interpolator;

		virtual std::pair< std::pair< double, double >, std::pair< double, double>>
		getInterpolationPoints(
				std::vector< std::pair< double, double > >::const_iterator & it ) const;
	};

	/// Class that behaves like suction curve. It is standard (linear or logarithmic) interpolation
	/// except for the results in first range where curve will return constant value equal to the
	/// first point
	class SuctionCurve : public Curve {
	public:
		SuctionCurve( const std::vector< std::pair< double, double > > & values,
									const Interpolator & interpolator = Interpolation::Logarithmic );

		SuctionCurve( const std::initializer_list< std::pair< double, double > > & list,
									const Interpolator & interpolator );

	protected:
		std::pair< std::pair< double, double >, std::pair< double, double > >
		getInterpolationPoints(
				std::vector< std::pair< double, double > >::const_iterator & it ) const override;

	};

	/// Simple constant curve.
	class Constant : public ICurve {
	public:
		Constant( const double value );

		double value( const double t_position = 0 ) const override;

	private:
		double m_Value;
	};

}