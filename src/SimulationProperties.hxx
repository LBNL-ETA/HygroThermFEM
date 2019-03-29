#pragma once

namespace HygroThermFEM
{
    class SimulationProperties
    {
    public:
        static SimulationProperties & Instance();

        void setIterationParameters(double t_RelaxationParameter,
                                    double t_ErrorTolerance,
                                    size_t t_MaxNumberOfIterations);

        double relaxationParamter() const;
        double errorTolerance() const;
        size_t maxNumberOfIterations() const;

    private:
        SimulationProperties();
        ~SimulationProperties() = default;

        double m_RelaxationParameter;
        double m_ErrorTolerance;
        size_t m_MaxNumberOfIterations;
    };
}   // namespace HygroThermFEM
