#include <cmath>
#include <iostream>
#include <utility>

#include "Material.hxx"
#include "State.hxx"
#include "Node2D.hxx"

namespace HygroThermFEM
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMaterial
    ///////////////////////////////////////////////////////////////////////////////////////////////
    IMaterial::IMaterial(std::string cs,
                         const double density,
                         const double porosity,
                         const double heatCapacity,
                         const double diffusionResistanceFactor,
                         const std::vector<std::pair<double, double>> & thermalConductivity,
                         const std::vector<std::pair<double, double>> & liquidTransportationCurve,
                         const std::vector<std::pair<double, double>> & sorptionCurve,
                         const double emissivity) :
        m_Name(std::move(cs)),
        m_Density(density),
        m_Porosity(porosity),
        m_HeatCapacity(heatCapacity),
        m_DiffusionResistanceFactor(diffusionResistanceFactor),
        m_ThermalConductivity(new TabularFunction(thermalConductivity, Variable::water)),
        m_LiquidTransportCoefficient(new LiquidTransportationCurve(liquidTransportationCurve)),
        m_SorptionCurve(new TabularFunction(sorptionCurve, Variable::humidity)),
        m_Emissivity(emissivity)
    {}

    double IMaterial::density() const
    {
        return m_Density;
    }

    double IMaterial::heatCapacity() const
    {
        return m_HeatCapacity;
    }

    double IMaterial::porosity() const
    {
        return m_Porosity;
    }

    double IMaterial::diffusionResistanceFactor() const
    {
        return m_DiffusionResistanceFactor;
    }

    const std::vector<std::pair<double, double>> & IMaterial::thermalConductivity() const
    {
        return m_ThermalConductivity->getCurve();
    }

    const std::vector<std::pair<double, double>> & IMaterial::liquidTransportationCurve() const
    {
        return m_LiquidTransportCoefficient->getCurve();
    }

    const std::vector<std::pair<double, double>> & IMaterial::sorptionCurve() const
    {
        return m_SorptionCurve->getCurve();
    }

    std::string IMaterial::name() const
    {
        return m_Name;
    }

    double IMaterial::emissivity() const
    {
        return m_Emissivity;
    }

    bool operator<(const IMaterial & lhs, const IMaterial & rhs)
    {
        return lhs.m_Name < rhs.m_Name;
    }

    bool operator>(const IMaterial & lhs, const IMaterial & rhs)
    {
        return rhs < lhs;
    }

    bool operator<=(const IMaterial & lhs, const IMaterial & rhs)
    {
        return !(rhs < lhs);
    }

    bool operator>=(const IMaterial & lhs, const IMaterial & rhs)
    {
        return !(lhs < rhs);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IGas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    IGas::IGas(const std::string & cs,
               double density,
               double porosity,
               double heatCapacity,
               double diffusionResistanceFactor,
               const std::vector<std::pair<double, double>> & thermalConductivity,
               const std::vector<std::pair<double, double>> & liquidTransportationCurve,
               const std::vector<std::pair<double, double>> & sorptionCurve,
               double emissivity,
               CavityStandard cavityStandard) :
        IMaterial(cs,
                  density,
                  porosity,
                  heatCapacity,
                  diffusionResistanceFactor,
                  thermalConductivity,
                  liquidTransportationCurve,
                  sorptionCurve,
                  emissivity),
        m_CavityStandard(cavityStandard)
    {}

    CavityStandard IGas::standard() const
    {
        return m_CavityStandard;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SolidMaterial
    ///////////////////////////////////////////////////////////////////////////////////////////////
    SolidMaterial::SolidMaterial(
      const std::string & name,
      const double density,
      const double porosity,
      const double heatCapacity,
      const double diffusionResistanceFactor,
      const std::vector<std::pair<double, double>> & thermalConductivity,
      const std::vector<std::pair<double, double>> & liquidTransportCurve,
      const std::vector<std::pair<double, double>> & sorptionCurve,
      const double emissivity) :
        IMaterial(name,
                  density,
                  porosity,
                  heatCapacity,
                  diffusionResistanceFactor,
                  thermalConductivity,
                  liquidTransportCurve,
                  sorptionCurve,
                  emissivity)
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

    double SolidMaterial::saturationConcentration(const INode2D & node)
    {
        const auto temperature = node.property(Variable::temperature);

        return saturationConcentrationAtTemperature(temperature);
    }

    double SolidMaterial::waterContent(const INode2D & node, WaterContent wContent) const
    {
        std::map<WaterContent, double> results;
        results[WaterContent::Water] = waterContent(node);
        results[WaterContent::Vapor] = vaporContent(node);
        results[WaterContent::Liquid] = liquidWaterContent(node);
        results[WaterContent::Ice] = iceContent(node);

        return results.at(wContent);
    }

    double SolidMaterial::waterContent(const INode2D & node) const
    {
        return m_SorptionCurve->value(node);
    }

    double SolidMaterial::vaporContent(const INode2D & node) const
    {
        return saturationConcentration(node) * airPorosity(node)
               * node.property(Variable::humidity);
    }

    double SolidMaterial::liquidWaterContent(const INode2D & node) const
    {
        return node.property(Variable::liquidPercent) * (waterContent(node) - vaporContent(node));
    }

    double SolidMaterial::iceContent(const INode2D & node) const
    {
        return (1 - node.property(Variable::liquidPercent))
               * (waterContent(node) - vaporContent(node));
    }

    double SolidMaterial::liquidPorosity(const INode2D & node) const
    {
        const auto waterContent = m_SorptionCurve->value(node);
        const auto maxWaterContent = m_SorptionCurve->maxY();
        return waterContent / maxWaterContent * m_Porosity;
    }

    double SolidMaterial::airPorosity(const INode2D & node) const
    {
        return m_Porosity - liquidPorosity(node);
    }

    void SolidMaterial::updateThermalConductivity(double)
    {
        // No thermal conductivity update for solid materials
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Gas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    double Gas::waterContent(const INode2D &, WaterContent) const
    {
        return 0;
    }

    void Gas::updateThermalConductivity(double thermalConductivity)
    {
        const auto minX = m_ThermalConductivity->minX();
        const auto maxX = m_ThermalConductivity->maxX();
        auto & curve = m_ThermalConductivity->getCurve();
        curve.clear();
        curve.emplace_back(minX, thermalConductivity);
        curve.emplace_back(maxX, thermalConductivity);
    }

    Gas::Gas(const std::string & name, CavityStandard cavityStandard) :
        IGas(name,
             0.0,              // Density
             1.0,              // Porosity
             0.0,              // Heat Capacity
             2,                // Diffusion resistance factor
             {{0.0, 1.8}},     // Thermal conductivity
             {{0.0, 0.0}},     // Liquid transportation curve
             {{0.0, 0.0}},     // Sorption curve
             0.0,              // emissivity
             cavityStandard)   // Standard used for cavity calculations

    {}
}   // namespace HygroThermFEM
