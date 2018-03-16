#include <cmath>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"

namespace MoisThermFEM {

	MultiDomain::MultiDomain() : m_ThermalDomain( Property::temperature ),
															 m_MoistureDomain( Property::humidity ) {}

	Solution
	MultiDomain::transient( std::vector< double > & temperature, std::vector< double > & humidity,
													const double t_DTime ) {
		auto temperatureError = std::numeric_limits< double >::max();
		auto humidityError = std::numeric_limits< double >::max();
		auto currentTemperature = temperature;
		auto currentHumidity = humidity;

		while ( temperatureError > ConvergenceError || humidityError > ConvergenceError ) {
			m_ThermalDomain.updateNodeValues( currentHumidity, Property::humidity );
			auto temperatureSolution = m_ThermalDomain.transient( temperature, t_DTime );
			temperatureError = normError( temperatureSolution, currentTemperature );
			currentTemperature = temperatureSolution;

			m_MoistureDomain.updateNodeValues( currentTemperature, Property::temperature );
			auto humiditySolution = m_MoistureDomain.transient( humidity, t_DTime );
			humidityError = normError( humiditySolution, currentHumidity );
			currentHumidity = humiditySolution;
		}

		return Solution{ currentTemperature, currentHumidity };
	}

	void MultiDomain::createElement( const Node2D & t_Node1, const Node2D & t_Node2,
																	 const Node2D & t_Node3, const Node2D & t_Node4,
																	 const Material & mat ) {

		m_ThermalDomain.elementsCreator().createThermalElement( t_Node1, t_Node2, t_Node3, t_Node4,
																														mat );
		m_MoistureDomain.elementsCreator().createMoistureElement( t_Node1, t_Node2, t_Node3, t_Node4,
																															mat );

	}

	void MultiDomain::createConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2,
																					const double t_ConvectionCoefficient, const double t_Porosity,
																					const double t_AirTemperature, const double t_Humidity ) {
		m_ThermalDomain.boundariesCreator().createConvectionBC( t_Node1, t_Node2,
																														t_ConvectionCoefficient,
																														t_AirTemperature );
		m_MoistureDomain.boundariesCreator().createMoistureBC( t_Node1, t_Node2,
																													 t_ConvectionCoefficient, t_Porosity,
																													 t_Humidity, t_AirTemperature );
	}

	void MultiDomain::createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp1,
																				 const double t_Temp2 ) {
		m_ThermalDomain.boundariesCreator().createTemperatureBC( t_Node1, t_Node2, t_Temp1, t_Temp2 );
	}

	void MultiDomain::createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp ) {
		m_ThermalDomain.boundariesCreator().createTemperatureBC( t_Node1, t_Node2, t_Temp );
	}

	void MultiDomain::createBlackBodyRadiationBC( const Node2D & t_Node1, const Node2D & t_Node2,
																								const double t_Emissivity,
																								const double t_RadiationTemperature ) {
		m_ThermalDomain.boundariesCreator().createBlackBodyRadiationBC( t_Node1, t_Node2, t_Emissivity,
																																		t_RadiationTemperature );
	}

	double
	MultiDomain::normError( const std::vector< double > & vec1, const std::vector< double > & vec2 ) {
		auto norm1 = norm( vec1 );
		auto norm2 = norm( vec2 );

		return std::abs( norm1 - norm2 ) / norm1;
	}

	Solution::Solution( const std::vector< double > & temperature,
											const std::vector< double > & humidity ) : temperature( temperature ),
																																 humidity( humidity ) {}
}