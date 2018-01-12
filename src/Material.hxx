#pragma once

#include <vector>

namespace MoisThermFEM {

	struct Material {
		Material() = delete;
		Material( double Density, double Porosity, double HeatCapacity, double ThermalConductivity,
							double DiffusionResistanceFactor, double LiquidTransportCoefficient,
							const std::vector< std::pair< double, double>> & SorptionCurve );

		double Density;
		double Porosity;
		double HeatCapacity;
		double ThermalConductivity;
		double DiffusionResistanceFactor;
		double LiquidTransportCoefficient;
		std::vector< std::pair< double, double > > SorptionCurve;
	};

}

