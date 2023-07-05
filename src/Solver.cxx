#include "Solver.hxx"
#include "Exceptions.hxx"
#include "SimulationProperties.hxx"
#include "NodePool.hxx"
#include "MultiDomain.hxx"
#include "FEMMath.hxx"
#include "FEMMatrices.hxx"
#include "LinearSolver.hxx"
#include "VectorOperators.hxx"
#include "ModelUtilities.hxx"
#include "TimestepData.hxx"

namespace HygroThermFEM
{
    Solution::Solution(const double dtime,
                       std::vector<double> temperature,
                       std::vector<double> humidity,
                       std::vector<double> waterContent,
                       std::vector<double> liquidWaterContent,
                       std::vector<double> vaporContent,
                       std::vector<double> iceContent,
                       std::vector<NodeFlux> heatFlux,
                       std::vector<NodeFlux> waterFlux,
                       const double temperatureError,
                       const double humidityError) :
        dTime(dtime),
        temperature(std::move(temperature)),
        humidity(std::move(humidity)),
        waterContent(std::move(waterContent)),
        liquidWaterContent(std::move(liquidWaterContent)),
        vaporContent(std::move(vaporContent)),
        iceContent(std::move(iceContent)),
        heatFlux(std::move(heatFlux)),
        waterFlux(std::move(waterFlux)),
        temperatureError(temperatureError),
        humidityError(humidityError)
    {}

    std::vector<double> HygroThermFEM::steadyState(SingleDomain & domain)
    {
        const auto B{steadyStateRightHandSide(domain)};
        const auto A{steadyStateLeftHandSide(domain)};
        return CLinearSolver::solveEigen(A, B);
    }

    Solution steadyState(HygroThermFEM::MultiDomain & domain)
    {
        const auto ConvergenceError = SimulationProperties::Instance().errorTolerance();
        auto temperatureError{std::numeric_limits<double>::max()};
        auto humidityError{std::numeric_limits<double>::max()};
        const auto MaxIterations = SimulationProperties::Instance().maxNumberOfIterations();
        size_t currentIteration{0};
        auto humidity = NodePool::Instance().properties(Variable::humidity);
        auto previousHumidity = humidity;
        auto temperature = NodePool::Instance().properties(Variable::temperature);
        auto previousTemperature = temperature;
        do
        {
            if(domain.simulateMoisture)
            {
                humidity = steadyState(domain.moistureDomain);
                humidityError = HygroThermFEM::errorNorm(humidity, previousHumidity);
                previousHumidity = humidity;
                NodePool::Instance().updateNodeValues(humidity, BaseVariable::humidity);
            }
            else
            {
                humidityError = 0;
            }
            if(domain.simulateThermal)
            {
                temperature = steadyState(domain.thermalDomain);
                temperatureError = HygroThermFEM::errorNorm(temperature, previousTemperature);
                previousTemperature = temperature;
                NodePool::Instance().updateNodeValues(temperature, BaseVariable::temperature);
            }
            else
            {
                temperatureError = 0;
            }
            ++currentIteration;
        } while(temperatureError > ConvergenceError || humidityError > ConvergenceError
                || currentIteration > MaxIterations);

        updateNodeValues(humidity, BaseVariable::humidity, true);
        updateNodeValues(temperature, BaseVariable::temperature, true);

        const auto waterContent = properties(Variable::water);
        const auto liquidContent = properties(Variable::liquid);
        const auto vaporContent = properties(Variable::vapor);
        const auto iceContent = properties(Variable::ice);

        const auto heatFlux = domain.thermalDomain.flux();
        const auto waterFlux = domain.moistureDomain.flux();

        return Solution{0,
                        temperature,
                        humidity,
                        waterContent,
                        liquidContent,
                        vaporContent,
                        iceContent,
                        heatFlux,
                        waterFlux,
                        temperatureError,
                        humidityError};
    }

    std::pair<std::vector<double>, bool>
      transientTimestep(SingleDomain & domain,
                        const std::vector<double> & previousTimestepStateValues,
                        double t_DTime,
                        size_t timestepIndex)
    {
        const auto RelaxParameter = SimulationProperties::Instance().relaxationParamter();
        const auto ConvergenceError = SimulationProperties::Instance().errorTolerance();
        const auto MaxIterations = SimulationProperties::Instance().maxNumberOfIterations();

        auto A = transientM_K_H_Matrix(domain, t_DTime, timestepIndex);
        auto B = transientMT_R_Vector(domain, previousTimestepStateValues, t_DTime, timestepIndex);

        std::vector<double> solution;
        bool converged{false};
        bool stopIterations{false};

        if(isLinear(domain))
        {
            solution = CLinearSolver::solveEigen(A, B);
            domain.postProcess(solution);
            converged = true;
        }
        else
        {
            solution = previousTimestepStateValues;
            std::vector<double> normSolution{previousTimestepStateValues};

            auto currentNorm = norm(solution);

            size_t numOfIterations = 0;

            while(!stopIterations && !converged)
            {
                const double previousNorm = currentNorm;
                auto temp = A * solution;
                temp = B - temp;

                auto dU = CLinearSolver::solveEigen(A, temp);

                solution = solution + dU * RelaxParameter;
                normSolution = solution + dU;

                domain.postProcess(solution);
                domain.postProcess(normSolution);

                currentNorm = norm(normSolution);

                ++numOfIterations;

                NodePool::Instance().updateNodeValues(
                  solution, domain.m_Property, domain.m_AutomaticUpdatePreviousTimestep);

                A = transientM_K_H_Matrix(domain, t_DTime, timestepIndex);
                B =
                  transientMT_R_Vector(domain, previousTimestepStateValues, t_DTime, timestepIndex);

                converged = (std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12))
                            <= ConvergenceError;

                stopIterations = numOfIterations > MaxIterations;
            }
        }

        updateNodeValues(solution, domain.m_Property, domain.m_AutomaticUpdatePreviousTimestep);

        return std::make_pair(solution, converged);
    }

    SingleTimestepSolution HygroThermFEM::transient(SingleDomain & domain,
                                                    const std::vector<double> & currentStateValues,
                                                    double t_DTime,
                                                    size_t timestepIndex)
    {
        std::vector<double> solution;
        bool converged{false};
        auto currentDivisionLevel{0u};
        auto maxDivisionLevel{Timesteps::Settings::Instance().getMaxDivisions()};
        double currentDTime{t_DTime};
        double totalTime{0};
        auto stateVariables{currentStateValues};
        unsigned numberOfSubtimesteps{Timesteps::Settings::Instance().getNumberOfSubtimesteps()};

        // In case program failed to converge, it will cut down step to smaller one and will perform
        // multiple consecutive simulations in order to achieve solution at requested timestep.
        while(totalTime < t_DTime)
        {
            domain.notify(currentDivisionLevel, unsigned(totalTime / currentDTime));
            std::tie(solution, converged) = transientTimestep(domain, stateVariables, currentDTime, timestepIndex);
            if(!converged)
            {
                currentDTime = currentDTime / numberOfSubtimesteps;
                ++currentDivisionLevel;
                if(currentDivisionLevel > maxDivisionLevel)
                {
                    throw SolutionFailedToConvergeException();
                }
            }
            else
            {
                stateVariables = solution;
                totalTime += currentDTime;
            }
        }

        return {solution, t_DTime};
    }

    namespace
    {
        std::tuple<SingleTimestepSolution, double, std::vector<double>>
          executeThermalSimulation(HygroThermFEM::MultiDomain & domain,
                                   const std::vector<double> & currentTemperature,
                                   const std::vector<double> & previousTimestepTemperature,
                                   const double dTime,
                                   size_t timestepIndex)
        {
            if(!domain.simulateThermal)
            {
                return std::make_tuple(
                  SingleTimestepSolution{currentTemperature, dTime}, 0.0, currentTemperature);
            }
            auto newTemperatureSolution =
              transient(domain.thermalDomain, previousTimestepTemperature, dTime, timestepIndex);
            auto newTemperatureError =
              HygroThermFEM::errorNorm(newTemperatureSolution.solution, currentTemperature);
            auto newCurrentTemperature = newTemperatureSolution.solution;
            return std::make_tuple(
              newTemperatureSolution, newTemperatureError, newCurrentTemperature);
        }

        std::tuple<SingleTimestepSolution, double, std::vector<double>>
          executeMoistureSimulation(HygroThermFEM::MultiDomain & domain,
                                    const std::vector<double> & currentHumidity,
                                    const std::vector<double> & previousTimestepHumidity,
                                    const double dTime,
                                    size_t timestepIndex)
        {
            if(!domain.simulateMoisture)
            {
                return std::make_tuple(
                  SingleTimestepSolution{currentHumidity, dTime}, 0.0, currentHumidity);
            }
            auto newHumiditySolution =
              transient(domain.moistureDomain, previousTimestepHumidity, dTime, timestepIndex);
            auto newHumidityError =
              HygroThermFEM::errorNorm(newHumiditySolution.solution, currentHumidity);
            auto newCurrentHumidity = newHumiditySolution.solution;
            return std::make_tuple(newHumiditySolution, newHumidityError, newCurrentHumidity);
        }

        void executeTransientIteration(HygroThermFEM::MultiDomain & domain,
                                       const std::vector<double> & previousTimestepTemperature,
                                       const std::vector<double> & previousTimestepHumidity,
                                       double & temperatureError,
                                       double & humidityError,
                                       std::vector<double> & currentTemperature,
                                       std::vector<double> & currentHumidity,
                                       SingleTimestepSolution & temperatureSolution,
                                       SingleTimestepSolution & humiditySolution,
                                       const double dTime,
                                       size_t timestepIndex,
                                       const double ConvergenceError,
                                       const size_t MaxIterations)
        {
            size_t localIterCounter{0};
            while(humidityError > ConvergenceError && temperatureError > ConvergenceError
                  && localIterCounter <= MaxIterations)
            {
                std::tie(humiditySolution, humidityError, currentHumidity) =
                  executeMoistureSimulation(
                    domain, currentHumidity, previousTimestepHumidity, dTime, timestepIndex);
                std::tie(temperatureSolution, temperatureError, currentTemperature) =
                  executeThermalSimulation(
                    domain, currentTemperature, previousTimestepTemperature, dTime, timestepIndex);
                ++localIterCounter;
            }
        }
    }   // namespace

    Solution transient(HygroThermFEM::MultiDomain & domain,
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

        // Note that temperature and humidity are solved separately first and then updated with new
        // data for next iteration.
        do
        {
            executeTransientIteration(domain,
                                      previousTimestepTemperature,
                                      previousTimestepHumidity,
                                      temperatureError,
                                      humidityError,
                                      currentTemperature,
                                      currentHumidity,
                                      temperatureSolution,
                                      humiditySolution,
                                      dTime,
                                      timestepIndex,
                                      ConvergenceError,
                                      MaxIterations);

            if(domain.simulateMoisture)
            {
                updateNodeValues(temperatureSolution.solution, BaseVariable::temperature, false);
                std::tie(humiditySolution, humidityError, currentHumidity) =
                  executeMoistureSimulation(
                    domain, currentHumidity, previousTimestepHumidity, dTime, timestepIndex);
            }

            if(domain.simulateThermal)
            {
                updateNodeValues(humiditySolution.solution, BaseVariable::humidity, false);
                std::tie(temperatureSolution, temperatureError, currentTemperature) =
                  executeThermalSimulation(
                    domain, currentTemperature, previousTimestepTemperature, dTime, timestepIndex);
            }

            ++currentIteration;
        }

        while((temperatureError > ConvergenceError && humidityError > ConvergenceError)
              || currentIteration > MaxIterations);

        updateNodeValues(temperatureSolution.solution, BaseVariable::temperature, true);
        updateNodeValues(humiditySolution.solution, BaseVariable::humidity, true);

        updateNodeValues(currentHumidity, BaseVariable::humidity, true);
        updateNodeValues(currentTemperature, BaseVariable::temperature, true);

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
}   // namespace HygroThermFEM