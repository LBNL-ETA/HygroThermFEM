#include "Material.hxx"
#include "Curve.hxx"
#include "FEMunique.hxx"

namespace MoisThermFEM {

	Material::Material( const std::string & Name, double Density, double Porosity, double HeatCapacity,
											double ThermalConductivity, double DiffusionResistanceFactor,
											const std::vector< std::pair< double, double > > & LiquidTransportCurve,
											const std::vector< std::pair< double, double > > & SorptionCurve ) :
			m_Name( Name ), m_Density( Density ), m_Porosity( Porosity ), m_HeatCapacity( HeatCapacity ),
			m_ThermalConductivity( ThermalConductivity ),
			m_DiffusionResistanceFactor( DiffusionResistanceFactor ),
			m_LiquidTransportCoefficient(
					fem::make_unique< FenestrationCommon::SuctionCurve >( LiquidTransportCurve ) ),
			m_SorptionCurve( fem::make_unique< FenestrationCommon::Curve >( SorptionCurve ) ) {}

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

	double Material::liquidTransportationCoefficient( const double t_value ) const {
		return m_LiquidTransportCoefficient->value( t_value );
	}

	double Material::sorption( const double t_value ) const {
		return m_SorptionCurve->firstDerivative( t_value );
	}

	std::string Material::name() const {
		return m_Name;
	}
}