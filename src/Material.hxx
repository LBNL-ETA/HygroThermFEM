#pragma once

#include <vector>
#include <memory>
#include <string>

#include "Functions.hxx"

namespace MoisThermFEM
{
	/// Used to simplify interface for certain classes. Instead of creating four different function,
	/// it is sufficient only to call one when passing this enum argument.
	enum class WaterContent{Water, Liquid, Vapor, Ice};

    class Material
    {
        friend class MaterialPool;

    public:
        Material() = delete;

        std::string name() const;
        double density() const;
        double heatCapacity() const;
        double porosity() const;
        double thermalConductivity() const;
        double diffusionResistanceFactor() const;
        std::vector<std::pair<double, double>> liquidTransportationCurve() const;

        std::vector<double> waterContent(const std::vector<double> & humidity) const;

        double waterContent(const State & t_State, WaterContent waterContent) const;

        std::vector<std::pair<double, double>> sorptionCurve() const;

        /// Materials will be stored in set which require operator >.
        friend bool operator<(const Material & lhs, const Material & rhs);
        friend bool operator>(const Material & lhs, const Material & rhs);
        friend bool operator<=(const Material & lhs, const Material & rhs);
        friend bool operator>=(const Material & lhs, const Material & rhs);

    private:
		double waterContent(const State & t_State) const;
		double vaporContent(const State & t_State) const;
		double liquidWaterContent(const State & t_State) const;
		double iceContent(const State & t_State) const;

        /// Create material by using MaterialPool.
        Material(const std::string & Name,
                 double Density,
                 double Porosity,
                 double HeatCapacity,
                 double ThermalConductivity,
                 double DiffusionResistanceFactor,
                 const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                 const std::vector<std::pair<double, double>> & SorptionCurve);

        std::string m_Name;
        double m_Density;
        double m_Porosity;
        double m_HeatCapacity;
        double m_ThermalConductivity;
        double m_DiffusionResistanceFactor;
        std::unique_ptr<MoisThermFEM::TabularFunction> m_LiquidTransportCoefficient;
        std::unique_ptr<MoisThermFEM::TabularFunction> m_SorptionCurve;

		double saturatedVaporContent(const State & t_State) const;
        double liquidPorosity(const State & t_State) const;
        double airPorosity(const State & t_State) const;
    };

}   // namespace MoisThermFEM
