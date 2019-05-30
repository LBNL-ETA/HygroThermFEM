#pragma once

#include <cstddef>

namespace HygroThermFEM
{
    class SimulationProperties
    {
    public:
        static SimulationProperties & Instance();

        void setIterationParameters(double t_RelaxationParameter,
                                    double t_ErrorTolerance,
                                    std::size_t t_MaxNumberOfIterations);

        double relaxationParamter() const;
        double errorTolerance() const;
        std::size_t maxNumberOfIterations() const;

    private:
        SimulationProperties();
        ~SimulationProperties() = default;

        struct DefaultProperties
        {
            double relaxationParameter{1};
            double errorTolerance{1e-5};
            double maxIterations{50};
        } defaultProperties;

        double m_RelaxationParameter;
        double m_ErrorTolerance;
        size_t m_MaxNumberOfIterations;
    };
}   // namespace HygroThermFEM
