#pragma once

#include <vector>
#include <memory>
#include <string>

#include "Functions.hxx"


namespace MoisThermFEM {

	class Material {

		friend class MaterialPool;

	public:
		Material() = delete;

		std::string name() const;

		double density() const;

		double heatCapacity() const;

		double porosity() const;

		double thermalConductivity() const;

		double diffusionResistanceFactor() const;

		std::vector< std::pair< double, double > > liquidTransportationCurve() const;

		std::vector< double > waterContent( const std::vector< double > & humidity ) const;
		double waterContent(const double humidity) const;
		std::vector< std::pair< double, double > > sorptionCurve() const;

	private:
		/// Create material by using MaterialPool.
		Material( const std::string & Name, double Density, double Porosity, double HeatCapacity,
							double ThermalConductivity,
							double DiffusionResistanceFactor,
							const std::vector< std::pair< double, double > > & LiquidTransportCurve,
							const std::vector< std::pair< double, double > > & SorptionCurve );

		std::string m_Name;
		double m_Density;
		double m_Porosity;
		double m_HeatCapacity;
		double m_ThermalConductivity;
		double m_DiffusionResistanceFactor;
		std::unique_ptr< MoisThermFEM::TabularFunction > m_LiquidTransportCoefficient;
		std::unique_ptr< MoisThermFEM::TabularFunction > m_SorptionCurve;
	};

}

