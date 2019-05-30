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
        m_MaxNumberOfIterations(defaultProperties.maxIterations)
    {}

    void SimulationProperties::setIterationParameters(const double t_RelaxationParameter,
                                                      const double t_ErrorTolerance,
                                                      const size_t t_MaxNumberOfIterations)
    {
        m_RelaxationParameter = t_RelaxationParameter;
        m_ErrorTolerance = t_ErrorTolerance;
        m_MaxNumberOfIterations = t_MaxNumberOfIterations;
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
