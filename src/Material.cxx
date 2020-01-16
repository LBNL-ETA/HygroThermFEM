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
      double thermalConductivityDry,
      double density,
      double porosity,
      double heatCapacity,
      double diffusionResistanceFactor,
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      double moistureDependentMeasurementTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      double temperatureDependentMeasurementHumidity,
      const std::vector<FenestrationCommon::point> & liquidTransportationCurve,
      const std::vector<FenestrationCommon::point> & sorptionCurve,
      double emissivity,
      const bool isLinear) :
        m_Name(std::move(name)),
        m_ThermalConductivityDry(thermalConductivityDry),
        m_Density(density),
        m_Porosity(porosity),
        m_SpecificHeatCapacity(heatCapacity),
        m_DiffusionResistanceFactor(diffusionResistanceFactor),
        m_ThermalConductivity2DTable(nullptr),
        m_LiquidTransportCoefficient(nullptr),
        m_SorptionCurve(nullptr),
        m_Emissivity(emissivity),
        m_Linear(isLinear)
    {
        setThermalConductivityMoistureAndTemperatureDependent(
          thermalConductivityMoistureDependent,
          moistureDependentMeasurementTemperature,
          thermalConductivityTemperatureDependent,
          temperatureDependentMeasurementHumidity);

        setLiquidTransportationCurve(liquidTransportationCurve);
        setSorptionCurve(sorptionCurve);
    }

    IMaterial::IMaterial(std::string name, const bool isLinear) :
        m_Name(std::move(name)), m_Linear(isLinear)
    {}

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

    void IMaterial::setThermalConductivity(const double thermalConductivity)
    {
        m_ThermalConductivityDry = thermalConductivity;
    }

    bool IMaterial::hasThermalConductivityDry() const
    {
        return m_ThermalConductivityDry.has_value();
    }

    double IMaterial::density() const
    {
        if(!m_Density.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have density.");
        }
        return m_Density.value();
    }

    void IMaterial::setDensity(const double density)
    {
        m_Density = density;
    }

    bool IMaterial::hasDensity() const
    {
        return m_Density.has_value();
    }

    double IMaterial::heatCapacity() const
    {
        if(!m_SpecificHeatCapacity.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have specific heat capacity.");
        }
        return m_SpecificHeatCapacity.value();
    }

    void IMaterial::setHeatCapacity(const double heatCapacity)
    {
        m_SpecificHeatCapacity = heatCapacity;
    }

    bool IMaterial::hasHeatCapacity() const
    {
        return m_SpecificHeatCapacity.has_value();
    }

    double IMaterial::porosity() const
    {
        if(!m_Porosity.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have porosity.");
        }
        return m_Porosity.value();
    }

    void IMaterial::setPorosity(const double porosity)
    {
        m_Porosity = porosity;
    }

    bool IMaterial::hasPorosity() const
    {
        return m_Porosity.has_value();
    }

    double IMaterial::emissivity() const
    {
        if(!m_Emissivity.has_value())
        {
            throw std::runtime_error("Material " + m_Name + " do not have emissivity.");
        }
        return m_Emissivity.value();
    }

    void IMaterial::setEmissivity(const double emissivity)
    {
        m_Emissivity = emissivity;
    }

    bool IMaterial::hasEmissivity() const
    {
        return m_Emissivity.has_value();
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

    void IMaterial::setDiffusionResistanceFactor(const double diffusionResistanceFactor)
    {
        m_DiffusionResistanceFactor = diffusionResistanceFactor;
    }

    bool IMaterial::hasDiffusionResistanceFactor() const
    {
        return m_DiffusionResistanceFactor.has_value();
    }

    bool IMaterial::isLinear() const
    {
        return m_Linear;
    }

    TabularFunction2D IMaterial::thermalConductivityMoistureAndTemperatureDependent() const
    {
        if(m_ThermalConductivity2DTable == nullptr)
        {
            throw std::runtime_error("Material " + m_Name
                                     + " do not have assigned table for moisture and temperature "
                                       "dependent thermal conductivity.");
        }
        return *m_ThermalConductivity2DTable;
    }

    void IMaterial::setThermalConductivityMoistureAndTemperatureDependent(
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      double moistureDependentMeasuredTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      double temperatureDependentMeasurementHumidity)
    {
        m_ThermalConductivity2DTable =
          std::make_unique<TabularFunction2D>(thermalConductivityMoistureDependent,
                                              moistureDependentMeasuredTemperature,
                                              Variable::humidity,
                                              thermalConductivityTemperatureDependent,
                                              temperatureDependentMeasurementHumidity,
                                              Variable::temperature);
    }

    bool IMaterial::hasThermalConductivityMoistureAndTemperatureDependent() const
    {
        return (m_ThermalConductivity2DTable != nullptr);
    }

    const std::vector<FenestrationCommon::point> & IMaterial::liquidTransportationCurve() const
    {
        if(m_LiquidTransportCoefficient == nullptr)
        {
            throw std::runtime_error("Material " + m_Name
                                     + " do not have liquid transportation coefficient.");
        }
        return m_LiquidTransportCoefficient->getCurve();
    }

    void IMaterial::setLiquidTransportationCurve(
      const std::vector<FenestrationCommon::point> & liquidTransportationCurve)
    {
        m_LiquidTransportCoefficient =
          std::make_unique<LiquidTransportationCurve>(liquidTransportationCurve);
    }

    bool IMaterial::hasLiquidTransportationCurve() const
    {
        return (m_LiquidTransportCoefficient != nullptr);
    }

    const std::vector<FenestrationCommon::point> & IMaterial::sorptionCurve() const
    {
        if(m_SorptionCurve == nullptr)
        {
            throw std::runtime_error("Material " + m_Name + " do not have sorption curve.");
        }
        return m_SorptionCurve->getCurve();
    }

    void IMaterial::setSorptionCurve(const std::vector<FenestrationCommon::point> & sorptionCurve)
    {
        m_SorptionCurve = std::make_unique<TabularFunction1D>(sorptionCurve, Variable::humidity);
    }

    bool IMaterial::hasSorptionCurve() const
    {
        return (m_SorptionCurve != nullptr);
    }

    double IMaterial::saturationConcentration(const INode2D & node)
    {
        const auto temperature = node.property(Variable::temperature);

        return saturationConcentrationAtTemperature(temperature);
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

    IGas::IGas(std::string name, const CavityStandard cavityStandard, Gases::CGas gas) :
        IMaterial(std::move(name), false), m_CavityStandard(cavityStandard), m_Gas(std::move(gas))
    {
        setThermalConductivity(DefaultThermalConductivityDry);
        setDensity(DefaultDensity);
        setPorosity(DefaultPorosity);
        setHeatCapacity(DefaultHeatCapacity);
        setDiffusionResistanceFactor(DefaultDiffusionResistanceFactor);
        setSorptionCurve(DefaultSorptionCurve);

        // Gas cavities in equations will use 2D table in order to estimate all necessary data.
        // That is why we need to insert some default in the first place.
        setThermalConductivityMoistureAndTemperatureDependent(
          DefaultThermalConductivityMoistureDependent,
          DefaultMoistureDependentMeasurementTemperature,
          DefaultThermalConductivityTemperatureDependent,
          DefaultTemperatureDependentMeasurementMoisture);
    }

    CavityStandard IGas::standard() const
    {
        return m_CavityStandard;
    }

    Gases::CGas & IGas::getGas()
    {
        return m_Gas;
    }

    const double IGas::DefaultThermalConductivityDry = 0.05;
    const double IGas::DefaultDensity = 0.0;
    const double IGas::DefaultPorosity = 1.0;
    const double IGas::DefaultHeatCapacity = 0.0;
    const double IGas::DefaultDiffusionResistanceFactor = 2;
    const std::vector<FenestrationCommon::point> IGas::DefaultSorptionCurve = {{0.0, 0.0},
                                                                               {1.0, 1.0}};
    const std::vector<FenestrationCommon::point> IGas::DefaultThermalConductivityMoistureDependent =
      {{0.0, 0.05}, {1.0, 0.05}};
    const double IGas::DefaultMoistureDependentMeasurementTemperature = 0;
    const std::vector<FenestrationCommon::point>
      IGas::DefaultThermalConductivityTemperatureDependent = {{0.0, 0.05}, {1.0, 0.05}};
    const double IGas::DefaultTemperatureDependentMeasurementMoisture = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // SolidMaterial
    ///////////////////////////////////////////////////////////////////////////////////////////////
    SolidMaterial::SolidMaterial(
      std::string name,
      double thermalConductivityDry,
      double density,
      double porosity,
      double heatCapacity,
      double diffusionResistanceFactor,
      const std::vector<FenestrationCommon::point> & thermalConductivityMoistureDependent,
      double moistureDependentMeasurementTemperature,
      const std::vector<FenestrationCommon::point> & thermalConductivityTemperatureDependent,
      double temperatureDependentMeasurementHumidity,
      const std::vector<FenestrationCommon::point> & liquidTransportCurve,
      const std::vector<FenestrationCommon::point> & sorptionCurve,
      const double emissivity) :
        IMaterial(std::move(name),
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
            if(m_ThermalConductivity2DTable->maxXFirstTable()
               != m_SorptionCurve
                    ->maxY())   // only moisture dependence must match to sorption curve
            {
                throw std::runtime_error("Thermal conductivity curve does not correspond to "
                                         "sorption curve. Maximum water "
                                         "content is not identical in both tables.");
            }

            if(m_LiquidTransportCoefficient->maxX() != m_SorptionCurve->maxY())
            {
                throw std::runtime_error(
                  "Liquid transportation coefficient table does not correspond to sorption "
                  "curve. "
                  "Maximum water content is not identical in both tables.");
            }
        }
        catch(const std::runtime_error & e)
        {
            std::cout << e.what();
        }
    }

    SolidMaterial::SolidMaterial(std::string name) : IMaterial(std::move(name))
    {}

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

    Water Gas::waterContent(const INode2D & node) const
    {
        const auto vaporConcentration = saturationConcentration(node)
               * node.property(Variable::humidity);
        // Gas have only vapor concentration
        return {vaporConcentration, 0, vaporConcentration, 0};
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

    void Gas::updateDiffusionResistanceFactor(const double diffusionResistanceFactor)
    {
        m_DiffusionResistanceFactor = diffusionResistanceFactor;
    }

    void Gas::updateSorptionCurve(double saturationTemperature)
    {
        auto & curve = m_SorptionCurve->getCurve();
        auto & last{curve[curve.size() - 1]};
        const auto vaporContent{
          saturationConcentrationAtTemperature(saturationTemperature - 273.15)};
        last.y = vaporContent;
    }

    Gas::Gas(const std::string & name, CavityStandard cavityStandard, Gases::CGas gas) : IGas(name, cavityStandard, gas)
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
