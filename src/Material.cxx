#include "Material.hxx"
#include "FEMunique.hxx"
#include "State.hxx"

namespace MoisThermFEM {

	Material::Material( const std::string & Name, double Density, double Porosity,
											double HeatCapacity,
											double ThermalConductivity, double DiffusionResistanceFactor,
											const std::vector< std::pair< double, double > > & LiquidTransportCurve,
											const std::vector< std::pair< double, double > > & SorptionCurve ) :
			m_Name( Name ), m_Density( Density ), m_Porosity( Porosity ), m_HeatCapacity( HeatCapacity ),
			m_ThermalConductivity( ThermalConductivity ),
			m_DiffusionResistanceFactor( DiffusionResistanceFactor ),
			m_LiquidTransportCoefficient(
					fem::make_unique< MoisThermFEM::SuctionFunction >( LiquidTransportCurve,
																														 Property::humidity ) ),
			m_SorptionCurve( nullptr ) {
		/// Need to add extra point at the end of sorption curve because iterations can go out of scope
		/// (can be greater than one) and that will produce incorrect results.
		std::vector< std::pair< double, double > > updatedResults{ SorptionCurve };

		/*const std::size_t size = updatedResults.size();
		auto x1 = updatedResults[ size - 2 ].first;
		auto x2 = updatedResults[ size - 1 ].first;
		auto y1 = updatedResults[ size - 2 ].second;
		auto y2 = updatedResults[ size - 1 ].second;

		auto der = ( y2 - y1 ) / ( x2 - x1 );
		double distance = 100;

		auto x3 = x2 + distance;
		auto y3 = y2 + der * distance;
		updatedResults.emplace_back( x3, y3 );*/

		m_SorptionCurve = fem::make_unique< MoisThermFEM::TabularFunction >( updatedResults,
																																				 Property::humidity );
	}

	double Material::density() const {
		return m_Density;
	}

	double Material::heatCapacity() const {
		return m_HeatCapacity;
	}

	double Material::porosity() const {
		return m_Porosity;
	}

	double Material::thermalConductivity() const {
		return m_ThermalConductivity;
	}

	double Material::diffusionResistanceFactor() const {
		return m_DiffusionResistanceFactor;
	}

	std::vector< std::pair< double, double > > Material::liquidTransportationCurve() const {
		return m_LiquidTransportCoefficient->getCurve();
	}

	std::vector< double > Material::waterContent( const std::vector< double > & humidity ) const {
		std::vector< double > result( humidity.size() );
		for ( auto i = 0u; i < humidity.size(); ++i ) {
			result[ i ] = m_SorptionCurve->value( State( 0, humidity[ i ], 0 ) );
		}
		return result;
	}

	std::vector< std::pair< double, double > > Material::sorptionCurve() const {
		return m_SorptionCurve->getCurve();
	}

	std::string Material::name() const {
		return m_Name;
	}
}