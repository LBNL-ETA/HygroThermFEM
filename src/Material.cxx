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
            result[i] = waterContent(humidity[i]);
        }
        return result;
    }

    double Material::waterContent(const double humidity) const
    {
        return m_SorptionCurve->value(State(0, humidity, 0));
    }

    std::vector<std::pair<double, double>> Material::sorptionCurve() const
    {
        return m_SorptionCurve->getCurve();
    }

    std::string Material::name() const
    {
        return m_Name;
    }
}   // namespace MoisThermFEM