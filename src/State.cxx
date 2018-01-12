#include "State.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////////////////////////
	///  State
	////////////////////////////////////////////////////////////////////////////

	State::State( const double t_Temperature, const double t_Humidity,
								const double t_Pressure ) {
		m_Property[ Property::temperature ] = t_Temperature;
		m_Property[ Property::humidity ] = t_Humidity;
		m_Property[ Property::pressure ] = t_Pressure;
	}

	double State::getValue( Property t_Property ) const {
		return m_Property.at( t_Property );
	}

	State::State( const State & other ) : m_Property( other.m_Property ) {

	}

	State & State::operator=( const State & other ) {
		m_Property = other.m_Property;
		return *this;
	}

	void State::setValue( Property t_Property, double t_Value ) {
		m_Property[ t_Property ] = t_Value;
	}

}