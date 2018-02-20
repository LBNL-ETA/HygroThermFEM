#include "State.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////////////////////////
	///  State
	////////////////////////////////////////////////////////////////////////////

	State::State( const double t_Temperature, const double t_Humidity,
								const double t_Pressure ) {
		m_Property[ Iteration::Current ][ Property::temperature ] = t_Temperature;
		m_Property[ Iteration::Current ][ Property::humidity ] = t_Humidity;
		m_Property[ Iteration::Current ][ Property::pressure ] = t_Pressure;

		m_Property[ Iteration::Previous ][ Property::temperature ] = t_Temperature;
		m_Property[ Iteration::Previous ][ Property::humidity ] = t_Humidity;
		m_Property[ Iteration::Previous ][ Property::pressure ] = t_Pressure;
	}

	double State::getValue( const Property t_Property, const Iteration t_Iteration ) const {
		return m_Property.at( t_Iteration ).at( t_Property );
	}

	void State::setValue( const Property t_Property, const double t_Value ) {
		m_Property[ Iteration::Previous ][ t_Property ] = m_Property[ Iteration::Current ][ t_Property ];
		m_Property[ Iteration::Current ][ t_Property ] = t_Value;
	}

	double State::getDeltaValue( const Property t_Property ) const {
		return m_Property.at( Iteration::Current ).at( t_Property ) -
					 m_Property.at( Iteration::Previous ).at( t_Property );
	}

	const std::map< Property, double > & State::getCurrentValues() const {
		return m_Property.at( Iteration::Current );
	}

}