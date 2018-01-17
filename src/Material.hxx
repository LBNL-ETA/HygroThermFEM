#pragma once

#include <vector>

namespace MoisThermFEM {

	struct Material {
		Material() = delete;
		Material( double Density, double Porosity, double HeatCapacity, double ThermalConductivity,
							double DiffusionResistanceFactor,
							const std::vector< std::pair< double, double > > & LiquidTransportCoefficient,
							const std::vector< std::pair< double, double > > & SorptionCurve );

		double Density;
		double Porosity;
		double HeatCapacity;
		double ThermalConductivity;
		double DiffusionResistanceFactor;
		std::vector< std::pair< double, double > > LiquidTransportCoefficient;
		std::vector< std::pair< double, double > > SorptionCurve;
	};

}

