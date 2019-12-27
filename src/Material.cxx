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
      std::string name,
      std::optional<double> thermalConductivityDry,
      std::optional<double> density,
      std::optional<double> porosity,
      std::optional<double> heatCapacity,
      std::optional<double> diffusionResistanceFactor,
      std::optional<std::vector<FenestrationCommon::point>> thermalConductivityMoistureDependent,
      std::optional<double> moistureDependentMeasurementTemperature,
      std::optional<std::vector<FenestrationCommon::point>> thermalConductivityTemperatureDependent,
      std::optional<double> temperatureDependentMeasurementHumidity,
      std::optional<std::vector<FenestrationCommon::point>> liquidTransportationCurve,
      std::optional<std::vector<FenestrationCommon::point>> sorptionCurve,
      std::optional<double> emissivity,
      const bool isLinear) :
        m_Name(std::move(name)),
        m_ThermalConductivityDry(std::move(thermalConductivityDry)),
        m_Density(std::move(density)),
        m_Porosity(std::move(porosity)),
        m_SpecificHeatCapacity(std::move(heatCapacity)),
        m_DiffusionResistanceFactor(std::move(diffusionResistanceFactor)),
        m_ThermalConductivity2DTable(nullptr),
        m_LiquidTransportCoefficient(nullptr),
        m_SorptionCurve(nullptr),
        m_Emissivity(std::move(emissivity)),
        m_Linear(isLinear)
    {
        if(thermalConductivityMoistureDependent.has_value()
           && moistureDependentMeasurementTemperature.has_value()
           && thermalConductivityTemperatureDependent.has_value()
           && temperatureDependentMeasurementHumidity.has_value())
        {
            m_ThermalConductivity2DTable =
              std::make_unique<TabularFunction2D>(thermalConductivityMoistureDependent.value(),
                                                  moistureDependentMeasurementTemperature.value(),
                                                  Variable::humidity,
                                                  thermalConductivityTemperatureDependent.value(),
                                                  temperatureDependentMeasurementHumidity.value(),
                                                  Variable::temperature);
        }
        if(liquidTransportationCurve.has_value())
        {
            m_LiquidTransportCoefficient =
              std::make_unique<LiquidTransportationCurve>(liquidTransportationCurve.value());
        }
        if(sorptionCurve.has_value())
        {
            m_SorptionCurve =
              std::make_unique<TabularFunction1D>(sorptionCurve.value(), Variable::humidity);
        }
    }

    std::string IMaterial::name() const
    {
        return m_Name;
    }

    double IMaterial::thermalConductivityDry() const
    {
        if(!m_ThermalConductivityDry.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have thermal conductivity.");
        }
        return m_ThermalConductivityDry.value();
    }

    double IMaterial::density() const
    {
        if(!m_Density.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have density.");
        }
        return m_Density.value();
    }

    double IMaterial::heatCapacity() const
    {
        if(!m_SpecificHeatCapacity.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have specific heat capacity.");
        }
        return m_SpecificHeatCapacity.value();
    }

    double IMaterial::porosity() const
    {
        if(!m_Porosity.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have porosity.");
        }
        return m_Porosity.value();
    }

    double IMaterial::emissivity() const
    {
        if(!m_Emissivity.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have emissivity.");
        }
        return m_Emissivity.value();
    }

    double IMaterial::diffusionResistanceFactor() const
    {
        if(!m_DiffusionResistanceFactor.has_value())
        {
            throw std::runtime_error("Material " + m_Name
                                     + " do not have diffusion resistance factor.");
        }
        return m_DiffusionResistanceFactor.value();
    }

    bool IMaterial::isLinear() const
    {
        return m_Linear;
    }

    TabularFunction2D IMaterial::thermalConductivityMoistureDependent() const
    {
        if(m_ThermalConductivity2DTable == nullptr)
        {
            throw std::runtime_error("Material " + m_Name
                                     + " do not have assigned table for moisture and temperature "
                                       "dependent thermal conductivity.");
        }
        return *m_ThermalConductivity2DTable;
    }

    const std::vector<FenestrationCommon::point> & IMaterial::liquidTransportationCurve() const
    {
        if(m_LiquidTransportCoefficient == nullptr)
        {
            throw std::runtime_error("Material " + m_Name + " do not have liquid transportation coefficient.");
        }
        return m_LiquidTransportCoefficient->getCurve();
    }

    const std::vector<FenestrationCommon::point> & IMaterial::sorptionCurve() const
    {
        if(m_SorptionCurve == nullptr)
        {
            throw std::runtime_error("Material " + m_Name + " do not have sorption curve.");
        }
        return m_SorptionCurve->getCurve();
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

    IGas::IGas(
      std::string name,
      std::optional<double> thermalConductivityDry,
      std::optional<double> density,
      std::optional<double> porosity,
      std::optional<double> heatCapacity,
      std::optional<double> diffusionResistanceFactor,
      std::optional<std::vector<FenestrationCommon::point>> thermalConductivityMoistureDependent,
      std::optional<double> moistureDependentMeasurementTemperature,
      std::optional<std::vector<FenestrationCommon::point>> thermalConductivityTemperatureDependent,
      std::optional<double> temperatureDependentMeasurementHumidity,
      std::optional<std::vector<FenestrationCommon::point>> liquidTransportationCurve,
      std::optional<std::vector<FenestrationCommon::point>> sorptionCurve,
      std::optional<double> emissivity,
      const CavityStandard cavityStandard) :
        IMaterial(std::move(name),
                  std::move(thermalConductivityDry),
                  std::move(density),
                  std::move(porosity),
                  std::move(heatCapacity),
                  std::move(diffusionResistanceFactor),
                  std::move(thermalConductivityMoistureDependent),
                  std::move(moistureDependentMeasurementTemperature),
                  std::move(thermalConductivityTemperatureDependent),
                  std::move(temperatureDependentMeasurementHumidity),
                  std::move(liquidTransportationCurve),
                  std::move(sorptionCurve),
                  std::move(emissivity),
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
      std::string name,
      std::optional<double> thermalConductivityDry,
      std::optional<double> density,
      std::optional<double> porosity,
      std::optional<double> heatCapacity,
      std::optional<double> diffusionResistanceFactor,
      std::optional<std::vector<FenestrationCommon::point>> thermalConductivityMoistureDependent,
      std::optional<double> moistureDependentMeasurementTemperature,
      std::optional<std::vector<FenestrationCommon::point>> thermalConductivityTemperatureDependent,
      std::optional<double> temperatureDependentMeasurementHumidity,
      std::optional<std::vector<FenestrationCommon::point>> liquidTransportCurve,
      std::optional<std::vector<FenestrationCommon::point>> sorptionCurve,
      const double emissivity) :
        IMaterial(std::move(name),
                  std::move(thermalConductivityDry),
                  std::move(density),
                  std::move(porosity),
                  std::move(heatCapacity),
                  std::move(diffusionResistanceFactor),
                  std::move(thermalConductivityMoistureDependent),
                  std::move(moistureDependentMeasurementTemperature),
                  std::move(thermalConductivityTemperatureDependent),
                  std::move(temperatureDependentMeasurementHumidity),
                  std::move(liquidTransportCurve),
                  std::move(sorptionCurve),
                  emissivity)
    {
        try
        {
            if(m_ThermalConductivity2DTable != nullptr && m_SorptionCurve != nullptr)
            {
                if(m_ThermalConductivity2DTable->maxXFirstTable()
                   != m_SorptionCurve
                        ->maxY())   // only moisture dependence must match to sorption curve
                {
                    throw std::runtime_error("Thermal conductivity curve does not correspond to "
                                             "sorption curve. Maximum water "
                                             "content is not identical in both tables.");
                }
            }

            if(m_LiquidTransportCoefficient != nullptr && m_SorptionCurve != nullptr)
            {
                if(m_LiquidTransportCoefficient->maxX() != m_SorptionCurve->maxY())
                {
                    throw std::runtime_error(
                      "Liquid transportation coefficient table does not correspond to sorption "
                      "curve. "
                      "Maximum water content is not identical in both tables.");
                }
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
        return waterContent / maxWaterContent * m_Porosity.value();
    }

    double SolidMaterial::airPorosity(const INode2D & node) const
    {
        return m_Porosity.value() - liquidPorosity(node);
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
             0.05,   // Thermal conductivity dry
             0.0,    // Density
             1.0,    // Porosity
             0.0,    // Heat Capacity
             2,      // Diffusion resistance factor
             std::optional<std::vector<FenestrationCommon::point>>(
               {{0.0, 0.05}, {1.0, 0.05}}),   // Thermal conductivity moisture dependent
             0,
             std::optional<std::vector<FenestrationCommon::point>>(
               {{0.0, 0.05}, {1.0, 0.05}}),   // Thermal conductivity temperature dependent
             0,
             std::optional<std::vector<FenestrationCommon::point>>(
               {{0.0, 0.0}, {1.0, 1.0}}),   // Liquid transportation curve
             std::optional<std::vector<FenestrationCommon::point>>(
               {{0.0, 0.0}, {1.0, 1.0}}),   // Sorption curve
             0.0,                           // emissivity
             cavityStandard)                // Standard used for cavity calculations

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
