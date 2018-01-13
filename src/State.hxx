#pragma once

#include <map>

namespace MoisThermFEM {

	enum class Property {
		temperature, humidity, pressure
	};

	// Simple class to hold state variables
	class State {
	public:
		explicit State( const double t_Temperature = 0, const double t_Humidity = 0,
										const double t_Pressure = 101325 );

		State( const State & other );

		State & operator=( const State & other );

		double getValue( const Property t_Property ) const;
		void setValue( const Property t_Property, const double t_Value );

	private:
		std::map< Property, double > m_Property;

	};

}