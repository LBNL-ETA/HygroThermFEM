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
	};

	/// Simple constant
	class Constant : public ICurve {
	public:
		Constant( const double value );

		double value( const double t_position = 0 ) const override;

	private:
		double m_Value;
	};

}