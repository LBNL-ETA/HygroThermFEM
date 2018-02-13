#pragma once

#include "BoundaryConditions2D.hxx"

namespace MoisThermFEM {

	class Node2D;

	class BoundaryCondition2DFactory : public BoundaryConditions2D {

	public:
		BoundaryCondition2DFactory() = default;

		void createConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
														 const double t_ConvectionCoefficient, const double t_AirTemperature );

		void createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp1,
															const double t_Temp2 );

		void createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp );

		void createBlackBodyRadiationBC( const Node2D & t_Node1, const Node2D & t_Node2,
																		 const double t_Emissivity, const double t_RadiationTemperature );

	};

}
