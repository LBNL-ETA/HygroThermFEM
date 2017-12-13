#pragma once

#include "Line2D.hxx"

namespace MoisThermFEM {

	struct Node2D;

	////////////////////////////////////////////////////////
	// ConvectionBC
	////////////////////////////////////////////////////////

	class ConvectionBC: public ILineLinear2D {
	public:
		ConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
		              const double t_ConvectionCoefficient, const double t_AirTemperature );

	protected:
		const double m_ConvectionCoefficient;
		const double m_AirTemperature;

	};

	////////////////////////////////////////////////////////
	// TemperatureBC
	////////////////////////////////////////////////////////

	// TemperatureBC will be just special case of convection BC with huge value
	// for film coefficients
	class TemperatureBC: public ConvectionBC {
	public:
		TemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_NodeTemperatures );
	};

}