#pragma once

#include "IBCLine2D.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////
	/// ConvectionBC
	////////////////////////////////////////////////////////

	class ConvectionBC : public IBCLinear2D {
	public:
		ConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
									const double t_ConvectionCoefficient, const double t_AirTemperature );

		virtual FenestrationCommon::Vector< double > R_Vector() const override;
		virtual FenestrationCommon::SquareMatrix< double > H_Matrix() const override;

	protected:
		const double m_ConvectionCoefficient;
		const double m_AirTemperature;

	};

	////////////////////////////////////////////////////////
	/// TemperatureBC
	////////////////////////////////////////////////////////

	// TemperatureBC will be just special case of convection BC with huge value
	// for film coefficients
	class TemperatureBC : public ConvectionBC {
	public:
		TemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_NodeTemperatures );
		TemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp1, const double t_Temp2 );
	};

	///////////////////////////////////////////////////////
	/// BlackBodyRadiationBC
	///////////////////////////////////////////////////////

	class BlackBodyRadiationBC : public IBCLinear2D {
	public:
		BlackBodyRadiationBC( const Node2D & t_Node1, const Node2D & t_Node2,
													const double t_Emissivity, const double t_RadiationTemperature );

		virtual FenestrationCommon::Vector< double > R_Vector() const override;
		virtual FenestrationCommon::SquareMatrix< double > H_Matrix() const override;

		/// DHMatrix seems unnecessary for now. Solution did converge without it.
		/// FenestrationCommon::SquareMatrix< double > D_HMatrix() const override;

	private:
		/// Radiative convective coefficient that needs to be calculated based on current temperatures
		FenestrationCommon::Vector< double > HRadiative() const;

		/// First derivative of radiative convection coefficient
		/// FenestrationCommon::Vector< double > DHRadiative() const;

		double m_RadiationTemperature;
		double m_Emissivity;
	};

}