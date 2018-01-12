#include "Material.hxx"

namespace MoisThermFEM {

	Material::Material( double Density, double Porosity, double HeatCapacity,
											double ThermalConductivity, double DiffusionResistanceFactor,
											double LiquidTransportCoefficient,
											const std::vector< std::pair< double, double > > & SorptionCurve ) :
			Density( Density ), Porosity( Porosity ), HeatCapacity( HeatCapacity ),
			ThermalConductivity( ThermalConductivity ),
			DiffusionResistanceFactor( DiffusionResistanceFactor ),
			LiquidTransportCoefficient( LiquidTransportCoefficient ), SorptionCurve( SorptionCurve ) {}
}