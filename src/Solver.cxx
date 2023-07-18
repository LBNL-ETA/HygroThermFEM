#include "Solver.hxx"

#include "NodePool.hxx"
#include "SingleDomain.hxx"
#include "MultiDomain.hxx"
#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    SingleDomainTransientSolver::SingleDomainTransientSolver(SingleDomain & domain) :
        m_Domain{domain}
    {}

    IterationResult SingleDomainTransientSolver::performDomainIteration(
      std::vector<double> & currentVariable,
      const std::vector<double> & previousTimestepVariable,
      double dTime,
      size_t timestepIndex)
    {
        {
            auto newValueSolution = transient(previousTimestepVariable, dTime, timestepIndex);
            auto newValueError =
              HygroThermFEM::errorNorm(newValueSolution.solution, currentVariable);
            updateNodeValues(newValueSolution.solution, baseVariableOf(m_Domain), false);
            currentVariable = newValueSolution.solution;
            return {newValueError, newValueSolution};
        }
    }

    TransientSolver::TransientSolver(MultiDomain & domain) : m_Domain{domain}
    {}

    Solution TransientSolver::transient(const std::vector<double> & previousTimestepTemperature,
                                        const std::vector<double> & previousTimestepHumidity,
                                        double t_DTime,
                                        size_t timestepIndex)
    {
        const auto ConvergenceError{SimulationProperties::Instance().errorTolerance()};
        auto temperatureError{std::numeric_limits<double>::max()};
        auto humidityError{std::numeric_limits<double>::max()};
        auto currentTemperature{previousTimestepTemperature};
        auto currentHumidity{previousTimestepHumidity};
        double dTime{t_DTime};
        SingleTimestepSolution temperatureSolution{previousTimestepTemperature, t_DTime};
        SingleTimestepSolution humiditySolution{previousTimestepHumidity, t_DTime};

        const auto MaxIterations = SimulationProperties::Instance().maxNumberOfIterations();

        size_t currentIteration{0};

        if(!m_TemperatureSolver && m_Domain.simulateThermal)
        {
            m_TemperatureSolver = createSolver(m_Domain.thermalDomain);
        }

        if(!m_HumiditySolver && m_Domain.simulateMoisture)
        {
            m_HumiditySolver = createSolver(m_Domain.moistureDomain);
        }

        do
        {
            if(m_Domain.simulateMoisture)
            {
                IterationResult moistureIteration = m_HumiditySolver->performDomainIteration(
                  currentHumidity, previousTimestepHumidity, dTime, timestepIndex);
                humidityError = moistureIteration.error;
                humiditySolution = moistureIteration.solution;
            }

            if(m_Domain.simulateThermal)
            {
                IterationResult thermalIteration = m_TemperatureSolver->performDomainIteration(
                  currentTemperature, previousTimestepTemperature, dTime, timestepIndex);
                temperatureError = thermalIteration.error;
                temperatureSolution = thermalIteration.solution;
            }

            ++currentIteration;
        } while((temperatureError > ConvergenceError && humidityError > ConvergenceError)
                || currentIteration > MaxIterations);

        updateNodeValues(temperatureSolution.solution, BaseVariable::temperature, true);
        updateNodeValues(humiditySolution.solution, BaseVariable::humidity, true);

        const auto waterContent{properties(Variable::water)};
        const auto liquidContent{properties(Variable::liquid)};
        const auto vaporContent{properties(Variable::vapor)};
        const auto iceContent{properties(Variable::ice)};

        const auto heatFlux = m_Domain.thermalDomain.flux();
        const auto waterFlux = m_Domain.moistureDomain.flux();

        return Solution{dTime,
                        currentTemperature,
                        currentHumidity,
                        waterContent,
                        liquidContent,
                        vaporContent,
                        iceContent,
                        heatFlux,
                        waterFlux,
                        temperatureError,
                        humidityError};
    }

    Solution TransientSolver::transient(const std::vector<double> & previousTimestepTemperature,
                                        const std::vector<double> & previousTimestepHumidity,
                                        double t_DTime)
    {
        return transient(previousTimestepTemperature, previousTimestepHumidity, t_DTime, 0);
    }
}   // namespace HygroThermFEM