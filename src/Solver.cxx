#include "Solver.hxx"

#include "NodePool.hxx"
#include "SingleDomain.hxx"
#include "MultiDomain.hxx"
#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    Solution TransientSolver::transient(MultiDomain & domain,
                                        const std::vector<double> & previousTimestepTemperature,
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

        do
        {
            if(domain.simulateMoisture)
            {
                IterationResult moistureIteration = performDomainIteration(domain.moistureDomain,
                                                                           currentHumidity,
                                                                           previousTimestepHumidity,
                                                                           dTime,
                                                                           timestepIndex);
                humidityError = moistureIteration.error;
                humiditySolution = moistureIteration.solution;
            }

            if(domain.simulateThermal)
            {
                IterationResult thermalIteration =
                  performDomainIteration(domain.thermalDomain,
                                         currentTemperature,
                                         previousTimestepTemperature,
                                         dTime,
                                         timestepIndex);
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

        const auto heatFlux = domain.thermalDomain.flux();
        const auto waterFlux = domain.moistureDomain.flux();

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

    Solution TransientSolver::transient(MultiDomain & domain,
                                        const std::vector<double> & previousTimestepTemperature,
                                        const std::vector<double> & previousTimestepHumidity,
                                        double t_DTime)
    {
        return transient(domain, previousTimestepTemperature, previousTimestepHumidity, t_DTime, 0);
    }

    IterationResult
      TransientSolver::performDomainIteration(SingleDomain & domain,
                                              std::vector<double> & currentVariable,
                                              const std::vector<double> & previousTimestepVariable,
                                              double dTime,
                                              size_t timestepIndex)
    {
        {
            auto newValueSolution = transient(domain, previousTimestepVariable, dTime, timestepIndex);
            auto newValueError =
              HygroThermFEM::errorNorm(newValueSolution.solution, currentVariable);
            updateNodeValues(newValueSolution.solution, baseVariableOf(domain), false);
            currentVariable = newValueSolution.solution;
            return {newValueError, newValueSolution};
        }
    }
}   // namespace HygroThermFEM