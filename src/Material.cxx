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
    IMaterial::IMaterial(
      std::string cs,
      const double thermalConductivityDry,
      const double density,
      const double porosity,
      const double heatCapacity,
      const double diffusionResistanceFactor,
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      const double moistureDependentMeasurementTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      const double temperatureDependentMeasurementHumidity,
      const std::vector<FenestrationCommon::point> & liquidTransportationCurve,
      const std::vector<FenestrationCommon::point> & sorptionCurve,
      const double emissivity,
      const bool isLinear) :
        m_Name(std::move(cs)),
        m_ThermalConductivityDry(thermalConductivityDry),
        m_Density(density),
        m_Porosity(porosity),
        m_HeatCapacity(heatCapacity),
        m_DiffusionResistanceFactor(diffusionResistanceFactor),
        m_LiquidTransportCoefficient(new LiquidTransportationCurve(liquidTransportationCurve)),
        m_SorptionCurve(new TabularFunction1D(sorptionCurve, Variable::humidity)),
        m_Emissivity(emissivity),
        m_Linear(isLinear)
    {
        m_ThermalConductivity2DTable =
          std::make_unique<TabularFunction2D>(thermalConductivityMoistureDependent,
                                              moistureDependentMeasurementTemperature,
                                              Variable::humidity,
                                              thermalConductivityTemperatureDependent,
                                              temperatureDependentMeasurementHumidity,
                                              Variable::temperature);
    }

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

    bool IMaterial::isLinear() const
    {
        return m_Linear;
    }

    const std::vector<FenestrationCommon::point> &
      IMaterial::thermalConductivityMoistureDependent() const
    {
        return m_ThermalConductivity2DTable->getCurve();
    }

    const std::vector<FenestrationCommon::point> & IMaterial::liquidTransportationCurve() const
    {
        return m_LiquidTransportCoefficient->getCurve();
    }

    const std::vector<FenestrationCommon::point> & IMaterial::sorptionCurve() const
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

    double IMaterial::thermalConductivityDry() const
    {
        return m_ThermalConductivityDry;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IGas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    IGas::IGas(
      const std::string & cs,
      const double thermalConductivityDry,
      const double density,
      const double porosity,
      const double heatCapacity,
      const double diffusionResistanceFactor,
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      const double moistureDependentMeasurementTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      const double temperatureDependentMeasurementHumidity,
      const std::vector<FenestrationCommon::point> & liquidTransportationCurve,
      const std::vector<FenestrationCommon::point> & sorptionCurve,
      const double emissivity,
      const CavityStandard cavityStandard) :
        IMaterial(cs,
                  thermalConductivityDry,
                  density,
                  porosity,
                  heatCapacity,
                  diffusionResistanceFactor,
                  thermalConductivityMoistureDependent,
                  moistureDependentMeasurementTemperature,
                  thermalConductivityTemperatureDependent,
                  temperatureDependentMeasurementHumidity,
                  liquidTransportationCurve,
                  sorptionCurve,
                  emissivity,
                  false),   // Gases introduce nonlinearty into domain
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
      const double thermalConductivityDry,
      const double density,
      const double porosity,
      const double heatCapacity,
      const double diffusionResistanceFactor,
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      double moistureDependentMeasurementTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      double temperatureDependentMeasurementHumidity,
      const std::vector<FenestrationCommon::point> & liquidTransportCurve,
      const std::vector<FenestrationCommon::point> & sorptionCurve,
      const double emissivity) :
        IMaterial(name,
                  thermalConductivityDry,
                  density,
                  porosity,
                  heatCapacity,
                  diffusionResistanceFactor,
                  thermalConductivityMoistureDependent,
                  moistureDependentMeasurementTemperature,
                  thermalConductivityTemperatureDependent,
                  temperatureDependentMeasurementHumidity,
                  liquidTransportCurve,
                  sorptionCurve,
                  emissivity)
    {
        try
        {
            if(m_ThermalConductivity2DTable->maxXFirstTable() != m_SorptionCurve->maxY()) // only moisture dependence must match to sorption curve
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

    Water SolidMaterial::waterContent(const INode2D & node) const
    {
        return {
          totalWaterContent(node), liquidWaterContent(node), vaporContent(node), iceContent(node)};
    }

    double SolidMaterial::totalWaterContent(const INode2D & node) const
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
        return node.property(Variable::liquidPercent)
               * (totalWaterContent(node) - vaporContent(node));
    }

    double SolidMaterial::iceContent(const INode2D & node) const
    {
        return (1 - node.property(Variable::liquidPercent))
               * (totalWaterContent(node) - vaporContent(node));
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

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Gas
    ///////////////////////////////////////////////////////////////////////////////////////////////

    Water Gas::waterContent(const INode2D &) const
    {
        return {};
    }

    void Gas::updateThermalConductivity(double thermalConductivity)
    {
        const auto minX = m_ThermalConductivity2DTable->minX();
        const auto maxX = m_ThermalConductivity2DTable->maxX();
        auto & curve = m_ThermalConductivity2DTable->getCurve();
        curve.clear();
        curve.emplace_back(minX, thermalConductivity);
        curve.emplace_back(maxX, thermalConductivity);
    }

    Gas::Gas(const std::string & name, CavityStandard cavityStandard) :
        IGas(name,
             0.05,                         // Thermal conductivity dry
             0.0,                          // Density
             1.0,                          // Porosity
             0.0,                          // Heat Capacity
             2,                            // Diffusion resistance factor
             {{0.0, 0.05}, {1.0, 0.05}},   // Thermal conductivity moisture dependent
             0,
             {{0.0, 0.05}, {1.0, 0.05}},   // Thermal conductivity temperature dependent
             0,
             {{0.0, 0.0}, {1.0, 1.0}},   // Liquid transportation curve
             {{0.0, 0.0}, {1.0, 1.0}},   // Sorption curve
             0.0,                        // emissivity
             cavityStandard)             // Standard used for cavity calculations

    {}

    Water::Water(double water, double liquid, double vapor, double ice) :
        m_Content{{WaterContent::Water, water},
                  {WaterContent::Liquid, liquid},
                  {WaterContent::Vapor, vapor},
                  {WaterContent::Ice, ice}}
    {}

    double Water::content(WaterContent content) const
    {
        return m_Content.at(content);
    }

    Water & Water::operator*(const double & other)
    {
        m_Content[WaterContent::Water] *= other;
        m_Content[WaterContent::Liquid] *= other;
        m_Content[WaterContent::Vapor] *= other;
        m_Content[WaterContent::Ice] *= other;
        return *this;
    }

    Water & Water::operator+=(const Water & other)
    {
        m_Content[WaterContent::Water] += other.content(WaterContent::Water);
        m_Content[WaterContent::Liquid] += other.content(WaterContent::Liquid);
        m_Content[WaterContent::Vapor] += other.content(WaterContent::Vapor);
        m_Content[WaterContent::Ice] += other.content(WaterContent::Ice);
        return *this;
    }

    Water & Water::operator/(const double & other)
    {
        m_Content[WaterContent::Water] /= other;
        m_Content[WaterContent::Liquid] /= other;
        m_Content[WaterContent::Vapor] /= other;
        m_Content[WaterContent::Ice] /= other;
        return *this;
    }
}   // namespace HygroThermFEM
