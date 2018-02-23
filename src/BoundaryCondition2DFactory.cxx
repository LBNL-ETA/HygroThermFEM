#include "BoundaryCondition2DFactory.hxx"
#include "BoundaryCondition2D.hxx"

namespace MoisThermFEM {

	void
	BoundaryCondition2DFactory::createConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
																									const double t_ConvectionCoefficient,
																									const double t_AirTemperature ) {
		m_BCs.push_back( fem::make_unique< ConvectionBC >( t_Node1, t_Node2, t_ConvectionCoefficient,
																											 t_AirTemperature ) );
	}

	void BoundaryCondition2DFactory::createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2,
																												const double t_Temp1,
																												const double t_Temp2 ) {
		m_BCs.push_back( fem::make_unique< TemperatureBC >( t_Node1, t_Node2, t_Temp1, t_Temp2 ) );
	}

	void BoundaryCondition2DFactory::createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2,
																												const double t_Temp ) {
		m_BCs.push_back( fem::make_unique< TemperatureBC >( t_Node1, t_Node2, t_Temp ) );
	}

	void BoundaryCondition2DFactory::createBlackBodyRadiationBC( const Node2D & t_Node1,
																															 const Node2D & t_Node2,
																															 const double t_Emissivity,
																															 const double t_RadiationTemperature ) {
		m_BCs.push_back( fem::make_unique< BlackBodyRadiationBC >( t_Node1, t_Node2, t_Emissivity,
																															 t_RadiationTemperature ) );
		m_Linear = false;
	}

	void BoundaryCondition2DFactory::createMoistureBC( const Node2D & t_Node1, const Node2D & t_Node2,
																										 const double t_ConvectiveCoefficient,
																										 const double t_AirHumidity ) {
		m_BCs.push_back( fem::make_unique< MoistureBC >( t_Node1, t_Node2, t_ConvectiveCoefficient,
																										 t_AirHumidity ) );
	}

	void BoundaryCondition2DFactory::createHumidityBC( const Node2D & t_Node1, const Node2D & t_Node2,
																										 const double t_Humidity ) {
		m_BCs.push_back( fem::make_unique< HumidityBC >( t_Node1, t_Node2, t_Humidity ) );
	}


}