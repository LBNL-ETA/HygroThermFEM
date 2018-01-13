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

	State::State( const State & other ) : m_Property( other.m_Property ) {

	}

	State & State::operator=( const State & other ) {
		m_Property = other.m_Property;
		return *this;
	}

	double State::getValue( const Property t_Property ) const {
		return m_Property.at( t_Property );
	}

	void State::setValue( const Property t_Property, const double t_Value ) {
		m_Property[ t_Property ] = t_Value;
	}

}