#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    SimulationProperties & HygroThermFEM::SimulationProperties::Instance()
    {
        static SimulationProperties m_Instance;
        return m_Instance;
    }

    bool SimulationProperties::excludeWaterLiquidTransportation() const
    {
        return m_ExcludeWaterLiquidTransportation;
    }

    bool SimulationProperties::excludeHeatOfEvaporation() const
    {
        return m_ExcludeHeatOfEvaporation;
    }

    bool SimulationProperties::excludeCapillaryConduction() const
    {
        return m_ExcludeCapillaryConduction;
    }

    bool SimulationProperties::excludeVaporDiffusionConduction() const
    {
        return m_ExcludeVaporDiffusionConduction;
    }

    void SimulationProperties::setIterationParameters(const double t_RelaxationParameter,
                                                      const double t_ErrorTolerance,
                                                      const size_t t_MaxNumberOfIterations)
    {
        m_RelaxationParameter = t_RelaxationParameter;
        m_ErrorToleranceMoisture = t_ErrorTolerance;
        m_MaxNumberOfIterations = t_MaxNumberOfIterations;
    }

    void SimulationProperties::resetIterationParameters()
    {
        m_RelaxationParameter = defaultProperties.relaxationParameter;
        m_ErrorToleranceMoisture = defaultProperties.errorToleranceMoisture;
        m_MaxNumberOfIterations = defaultProperties.maxIterations;
    }

    void SimulationProperties::setCalculationParameters(const bool excludeWaterLiquidTransportation,
                                                        const bool excludeHeatOfEvaporation,
                                                        const bool excludeCapillaryConduction,
                                                        const bool excludeVaporDiffusionConduction,
                                                        const bool thermalConductivityMoistureAndTemperatureDependent)
    {
        m_ExcludeWaterLiquidTransportation = excludeWaterLiquidTransportation;
        m_ExcludeHeatOfEvaporation = excludeHeatOfEvaporation;
        m_ExcludeCapillaryConduction = excludeCapillaryConduction;
        m_ExcludeVaporDiffusionConduction = excludeVaporDiffusionConduction;
        m_ThermalConductivityMoistureAndTemperatureDependent = thermalConductivityMoistureAndTemperatureDependent;
    }

    void SimulationProperties::resetCalculationParameters()
    {
        m_ExcludeWaterLiquidTransportation = defaultProperties.excludeWaterLiquidTransportation;
        m_ExcludeHeatOfEvaporation = defaultProperties.excludeHeatOfEvaporation;
        m_ExcludeCapillaryConduction = defaultProperties.excludeCapillaryConduction;
        m_ExcludeVaporDiffusionConduction = defaultProperties.excludeVaporDiffusionConduction;
        m_ThermalConductivityMoistureAndTemperatureDependent = defaultProperties.thermalConductivityMoistureAndTemperatureDependent;
    }

    void SimulationProperties::reset()
    {
        resetIterationParameters();
        resetCalculationParameters();
    }

    double SimulationProperties::relaxationParamter() const
    {
        return m_RelaxationParameter;
    }

    double SimulationProperties::errorToleranceMoisture() const
    {
        return m_ErrorToleranceMoisture;
    }

    double SimulationProperties::errorToleranceTemperature() const
    {
        return m_ErrorToleranceTemperature;
    }

    size_t SimulationProperties::maxNumberOfIterations() const
    {
        return m_MaxNumberOfIterations;
    }

    bool SimulationProperties::thermalConductivityTemperatureAndMoistureDependent() const {
        return m_ThermalConductivityMoistureAndTemperatureDependent;
    }

}   // namespace HygroThermFEM
