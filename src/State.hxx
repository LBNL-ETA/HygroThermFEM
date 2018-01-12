#pragma once

#include <map>

namespace MoisThermFEM {

	////////////////////////////////////////////////////////////////////////////
	///  State
	////////////////////////////////////////////////////////////////////////////
	enum class Property {
		temperature, humidity, pressure
	};

	class State {
	public:
		explicit State( const double t_Temperature = 0, const double t_Humidity = 0,
										const double t_Pressure = 101325 );

		State( const State & other );

		State & operator=( const State & other );

		double getValue( Property t_Property ) const;
		void setValue( Property t_Property, double t_Value );

	private:
		std::map< Property, double > m_Property;

	};

}