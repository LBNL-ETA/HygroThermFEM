#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    SimulationProperties & HygroThermFEM::SimulationProperties::Instance()
    {
        static SimulationProperties m_Instance;
        return m_Instance;
    }

    SimulationProperties::SimulationProperties() :
        m_RelaxationParameter(defaultProperties.relaxationParameter),
        m_ErrorTolerance(defaultProperties.errorTolerance),
        m_MaxNumberOfIterations(defaultProperties.maxIterations),
        m_ExcludeWaterLiquidTransportation(defaultProperties.excludeWaterLiquidTransportation),
        m_ExcludeHeatOfEvaporation(defaultProperties.excludeHeatOfEvaporation),
        m_ExcludeCapillaryConduction(defaultProperties.excludeCapillaryConduction),
        m_ExcludeVaporDiffusionConduction(defaultProperties.excludeVaporDiffusionConduction)
    {}

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
        m_ErrorTolerance = t_ErrorTolerance;
        m_MaxNumberOfIterations = t_MaxNumberOfIterations;
    }

    void SimulationProperties::resetIterationParameters()
    {
        m_RelaxationParameter = defaultProperties.relaxationParameter;
        m_ErrorTolerance = defaultProperties.errorTolerance;
        m_MaxNumberOfIterations = defaultProperties.maxIterations;
    }

    void SimulationProperties::setCalculationParameters(bool excludeWaterLiquidTransportation,
                                                        bool excludeHeatOfEvaporation,
                                                        bool excludeCapillaryConduction,
                                                        bool excludeVaporDiffusionConduction)
    {
        m_ExcludeWaterLiquidTransportation = excludeWaterLiquidTransportation;
        m_ExcludeHeatOfEvaporation = excludeHeatOfEvaporation;
        m_ExcludeCapillaryConduction = excludeCapillaryConduction;
        m_ExcludeVaporDiffusionConduction = excludeVaporDiffusionConduction;
    }

    void SimulationProperties::resetCalculationParameters()
    {
        m_ExcludeWaterLiquidTransportation = defaultProperties.excludeWaterLiquidTransportation;
        m_ExcludeHeatOfEvaporation = defaultProperties.excludeHeatOfEvaporation;
        m_ExcludeCapillaryConduction = defaultProperties.excludeCapillaryConduction;
        m_ExcludeVaporDiffusionConduction = defaultProperties.excludeVaporDiffusionConduction;
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

    double SimulationProperties::errorTolerance() const
    {
        return m_ErrorTolerance;
    }

    size_t SimulationProperties::maxNumberOfIterations() const
    {
        return m_MaxNumberOfIterations;
    }

}   // namespace HygroThermFEM
