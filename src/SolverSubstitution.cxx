#include "SolverSubstitution.hxx"
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

namespace HygroThermFEM::Substitution
{
    namespace
    {
        /// Post processing for iterative solver. This needs to be connected with a specific solver
        /// and it is not necessary to keep it for other solvers.

        using PostProcessFunc = std::function<void(std::vector<double> &)>;

        // Define a map from SingleDomainType to PostProcessFunc.
        std::map<DomainType, PostProcessFunc> postProcessFuncMap = {
          {DomainType::Thermal,
           [](std::vector<double> & solution) {
               for(auto & val : solution)
               {
                   if(val < Constants::ABSOLUTEZERO)
                   {
                       val = Constants::ABSOLUTEZERO + 1e-6;
                   }
                   if(val > 1000)
                   {
                       val = 1000;
                   }
               }
           }},
          {DomainType::Moisture, [](std::vector<double> & solution) {
               for(auto & val : solution)
               {
                   if(val > 1)
                   {
                       val = 1;
                   }
                   if(val < 0)
                   {
                       val = 0;
                   }
               }
           }}};

        void postProcess(SingleDomain & domain, std::vector<double> & solution)
        {
            if(!domain.gasCavities.has_value())
            {
                domain.gasCavities.emplace(domain.elements);
                domain.gasCavities->setGravityVector(domain.gravityVector);
            }
            domain.gasCavities->update();

            // Domain-specific processing.
            postProcessFuncMap[domain.domainType](solution);
        }
    }   // namespace


    std::vector<double> steadyState(SingleDomain & domain)
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
                HygroThermFEM::updateNodeValues(humidity, BaseVariable::humidity, false);
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
                HygroThermFEM::updateNodeValues(temperature, BaseVariable::temperature, false);
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

    namespace
    {
        //! \brief Calling timestep calculations
        //! @param previousTimestepStateValues Current state values from previous timestep
        //! @param t_DTime Time different for between timesteps
        //! @param timestepIndex Current timestep index used in variable boundary conditions
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
            auto B =
              transientMT_R_Vector(domain, previousTimestepStateValues, t_DTime, timestepIndex);

            std::vector<double> solution;
            bool converged{false};
            bool stopIterations{false};

            if(isLinear(domain))
            {
                solution = CLinearSolver::solveEigen(A, B);
                postProcess(domain, solution);
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

                    postProcess(domain, solution);
                    postProcess(domain, normSolution);

                    currentNorm = norm(normSolution);

                    ++numOfIterations;

                    updateNodeValues(solution, baseVariable(domain), true);

                    A = transientM_K_H_Matrix(domain, t_DTime, timestepIndex);
                    B = transientMT_R_Vector(
                      domain, previousTimestepStateValues, t_DTime, timestepIndex);

                    converged = (std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12))
                                <= ConvergenceError;

                    stopIterations = numOfIterations > MaxIterations;
                }
            }

            updateNodeValues(solution, baseVariable(domain), true);

            return std::make_pair(solution, converged);
        }
    }   // namespace

    SingleTimestepSolution transient(SingleDomain & domain,
                                     const std::vector<double> & previousTimestepValues,
                                     double t_DTime,
                                     size_t timestepIndex)
    {
        std::vector<double> solution;
        bool converged{false};
        auto currentDivisionLevel{0u};
        auto maxDivisionLevel{Timesteps::Settings::Instance().getMaxDivisions()};
        double currentDTime{t_DTime};
        double totalTime{0};
        auto stateVariables{previousTimestepValues};
        unsigned numberOfSubtimesteps{Timesteps::Settings::Instance().getNumberOfSubtimesteps()};

        // In case program failed to converge, it will cut down step to smaller one and will perform
        // multiple consecutive simulations in order to achieve solution at requested timestep.
        while(totalTime < t_DTime)
        {
            domain.notify(currentDivisionLevel, unsigned(totalTime / currentDTime));
            std::tie(solution, converged) =
              transientTimestep(domain, stateVariables, currentDTime, timestepIndex);
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
          executeSimulation(HygroThermFEM::SingleDomain & domain,
                            const std::vector<double> & currentValue,
                            const std::vector<double> & previousTimestepValue,
                            const double dTime,
                            size_t timestepIndex)
        {
            auto newValueSolution = transient(domain, previousTimestepValue, dTime, timestepIndex);
            auto newValueError = HygroThermFEM::errorNorm(newValueSolution.solution, currentValue);
            auto newCurrentValue = newValueSolution.solution;
            return std::make_tuple(newValueSolution, newValueError, newCurrentValue);
        }


        void executeTransientIteration(HygroThermFEM::MultiDomain & multiDomain,
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
                if(multiDomain.simulateMoisture)
                {
                    std::tie(humiditySolution, humidityError, currentHumidity) =
                      executeSimulation(multiDomain.moistureDomain,
                                        currentHumidity,
                                        previousTimestepHumidity,
                                        dTime,
                                        timestepIndex);
                }
                if(multiDomain.simulateThermal)
                {
                    std::tie(temperatureSolution, temperatureError, currentTemperature) =
                      executeSimulation(multiDomain.thermalDomain,
                                        currentTemperature,
                                        previousTimestepTemperature,
                                        dTime,
                                        timestepIndex);
                }
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
                  executeSimulation(domain.moistureDomain,
                                    currentHumidity,
                                    previousTimestepHumidity,
                                    dTime,
                                    timestepIndex);
            }

            if(domain.simulateThermal)
            {
                updateNodeValues(humiditySolution.solution, BaseVariable::humidity, false);
                std::tie(temperatureSolution, temperatureError, currentTemperature) =
                  executeSimulation(domain.thermalDomain,
                                    currentTemperature,
                                    previousTimestepTemperature,
                                    dTime,
                                    timestepIndex);
            }

            ++currentIteration;
        }

        while((temperatureError > ConvergenceError && humidityError > ConvergenceError)
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
}   // namespace HygroThermFEM