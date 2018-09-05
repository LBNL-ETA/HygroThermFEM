#include <cmath>

#include "Material.hxx"
#include "FEMunique.hxx"
#include "State.hxx"

namespace MoisThermFEM
{
    Material::Material(const std::string & Name,
                       double Density,
                       double Porosity,
                       double HeatCapacity,
                       double ThermalConductivity,
                       double DiffusionResistanceFactor,
                       const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                       const std::vector<std::pair<double, double>> & SorptionCurve) :
        m_Name(Name),
        m_Density(Density),
        m_Porosity(Porosity),
        m_HeatCapacity(HeatCapacity),
        m_ThermalConductivity(ThermalConductivity),
        m_DiffusionResistanceFactor(DiffusionResistanceFactor),
        m_LiquidTransportCoefficient(
          SuctionFunction::create(LiquidTransportCurve, Property::humidity)),
        m_SorptionCurve(TabularFunction::create(SorptionCurve, Property::humidity))
    {}

    bool operator<(const Material & lhs, const Material & rhs)
    {
        return lhs.m_Name < rhs.m_Name;
    }

    bool operator>(const Material & lhs, const Material & rhs)
    {
        return rhs < lhs;
    }

    bool operator<=(const Material & lhs, const Material & rhs)
    {
        return !(rhs < lhs);
    }

    bool operator>=(const Material & lhs, const Material & rhs)
    {
        return !(lhs < rhs);
    }

    double Material::density() const
    {
        return m_Density;
    }

    double Material::heatCapacity() const
    {
        return m_HeatCapacity;
    }

    double Material::porosity() const
    {
        return m_Porosity;
    }

    double Material::thermalConductivity() const
    {
        return m_ThermalConductivity;
    }

    double Material::diffusionResistanceFactor() const
    {
        return m_DiffusionResistanceFactor;
    }

    std::vector<std::pair<double, double>> Material::liquidTransportationCurve() const
    {
        return m_LiquidTransportCoefficient->getCurve();
    }

    std::vector<double> Material::waterContent(const std::vector<double> & humidity) const
    {
        std::vector<double> result(humidity.size());
        for(auto i = 0u; i < humidity.size(); ++i)
        {
            result[i] = waterContent(State(0, humidity[i], 0, 0));
        }
        return result;
    }

    double Material::saturatedVaporContent(const State & t_State) const
    {
        const auto temperature = t_State.getValue(Property::temperature);

        auto temp = 77.345 + 0.0057 * temperature - 7235.0 / temperature;
        temp = std::exp(temp);
        temp = temp / (461.4 * std::pow(temperature, 9.2));
        return temp;
    }

    double Material::waterContent(const State & t_State) const
    {
        return m_SorptionCurve->value(t_State);
    }

    double Material::vaporContent(const State & t_State) const
    {
        return saturatedVaporContent(t_State) * airPorosity(t_State)
               * t_State.getValue(Property::humidity);
    }

    double Material::liquidWaterContent(const State & t_State) const
    {
        return t_State.getLiquidPercent() * (waterContent(t_State) - vaporContent(t_State));
    }

    double Material::iceContent(const State & t_State) const
    {
        return (1 - t_State.getLiquidPercent()) * (waterContent(t_State) - vaporContent(t_State));
    }

    std::vector<std::pair<double, double>> Material::sorptionCurve() const
    {
        return m_SorptionCurve->getCurve();
    }

    std::string Material::name() const
    {
        return m_Name;
    }

    double Material::liquidPorosity(const State & t_State) const
    {
        const auto waterContent = m_SorptionCurve->value(t_State);
        const auto maxWaterContent = m_SorptionCurve->max();
        return waterContent / maxWaterContent * m_Porosity;
    }

    double Material::airPorosity(const State & t_State) const
    {
        return m_Porosity - liquidPorosity(t_State);
    }
}   // namespace MoisThermFEM