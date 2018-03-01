#pragma once

#include "Domains.hxx"

namespace MoisThermFEM {

	struct Solution {
		Solution( const std::vector< double > & temperature, const std::vector< double > & humidity );

		std::vector< double > temperature;
		std::vector< double > humidity;
	};

	class MultiDomain {
	public:
		MultiDomain();

		/// Calculates next timestep value from current values
		Solution
		transient( std::vector< double > & temperature, std::vector< double > & humidity, const double t_DTime );

		void createElement( const Node2D & t_Node1, const Node2D & t_Node2, const Node2D & t_Node3,
												const Node2D & t_Node4, const Material & mat );

		void createConvectionBC( const Node2D & t_Node1, const Node2D & t_Node2, const Material & material,
														 const double t_ConvectionCoefficient, const double t_AirTemperature,
														 const double t_Humidity );

		void createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp1,
															const double t_Temp2 );

		void createTemperatureBC( Node2D & t_Node1, Node2D & t_Node2, const double t_Temp );

		void createBlackBodyRadiationBC( const Node2D & t_Node1, const Node2D & t_Node2,
																		 const double t_Emissivity,
																		 const double t_RadiationTemperature );

	private:
		static double normError( const std::vector< double > & vec1, const std::vector< double > & vec2 );
		Domain m_ThermalDomain;
		Domain m_MoistureDomain;

	};

}
