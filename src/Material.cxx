#include <cmath>
#include <iostream>

#include "Material.hxx"
#include "State.hxx"
#include "Node2D.hxx"

namespace HygroThermFEM
{
    Material::Material(const std::string & Name,
                       const double Density,
                       const double Porosity,
                       const double HeatCapacity,
                       const double DiffusionResistanceFactor,
                       const std::vector<std::pair<double, double>> & ThermalConductivity,
                       const std::vector<std::pair<double, double>> & LiquidTransportCurve,
                       const std::vector<std::pair<double, double>> & SorptionCurve) :
        m_Name(Name),
        m_Density(Density),
        m_Porosity(Porosity),
        m_HeatCapacity(HeatCapacity),
        m_DiffusionResistanceFactor(DiffusionResistanceFactor),
        m_ThermalConductivity(new TabularFunction(ThermalConductivity, Variable::water)),
        m_LiquidTransportCoefficient(new SuctionCurve(LiquidTransportCurve)),
        m_SorptionCurve(new TabularFunction(SorptionCurve, Variable::humidity))
    {
        try
        {
            if(m_ThermalConductivity->maxX() != m_SorptionCurve->maxY())
            {
                throw std::runtime_error(
                  "Thermal conductivity curve does not correspond to sorption curve. Maximum water "
                  "content is not identical in both tables.");
            }
            if(m_LiquidTransportCoefficient->maxX() != m_SorptionCurve->maxY())
            {
                throw std::runtime_error(
                  "Liquid transportation coefficient table does not correspond to sorption curve. "
                  "Maximum water content is not identical in both tables.");
            }
        }
        catch(const std::runtime_error & e)
        {
            std::cout << e.what();
        }
    }

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

    double Material::diffusionResistanceFactor() const
    {
        return m_DiffusionResistanceFactor;
    }

    const std::vector<std::pair<double, double>> & Material::thermalConductivity() const
    {
        return m_ThermalConductivity->getCurve();
    }

    const std::vector<std::pair<double, double>> & Material::liquidTransportationCurve() const
    {
        return m_LiquidTransportCoefficient->getCurve();
    }

    double Material::saturatedVaporContent(const Node2D & node) const
    {
        const auto temperature = node.property(Variable::temperature);

        const auto tempinK = temperature + 273.15;

        auto temp = 77.345 + 0.0057 * tempinK - 7235.0 / tempinK;
        temp = std::exp(temp);
        temp = temp / (461.4 * std::pow(tempinK, 9.2));
        return temp;
    }

    double Material::waterContent(const Node2D & node, WaterContent wContent) const
    {
        std::map<WaterContent, double> results;
        results[WaterContent::Water] = waterContent(node);
        results[WaterContent::Vapor] = vaporContent(node);
        results[WaterContent::Liquid] = liquidWaterContent(node);
        results[WaterContent::Ice] = iceContent(node);

        return results.at(wContent);
    }

    double Material::waterContent(const Node2D & node) const
    {
        return m_SorptionCurve->value(node);
    }

    double Material::vaporContent(const Node2D & node) const
    {
        return saturatedVaporContent(node) * airPorosity(node) * node.property(Variable::humidity);
    }

    double Material::liquidWaterContent(const Node2D & node) const
    {
        return node.property(Variable::liquidPercent) * (waterContent(node) - vaporContent(node));
    }

    double Material::iceContent(const Node2D & node) const
    {
        return (1 - node.property(Variable::liquidPercent))
               * (waterContent(node) - vaporContent(node));
    }

    const std::vector<std::pair<double, double>> & Material::sorptionCurve() const
    {
        return m_SorptionCurve->getCurve();
    }

    std::string Material::name() const
    {
        return m_Name;
    }

    double Material::liquidPorosity(const Node2D & node) const
    {
        const auto waterContent = m_SorptionCurve->value(node);
        const auto maxWaterContent = m_SorptionCurve->maxY();
        return waterContent / maxWaterContent * m_Porosity;
    }

    double Material::airPorosity(const Node2D & node) const
    {
        return m_Porosity - liquidPorosity(node);
    }
}   // namespace HygroThermFEM
