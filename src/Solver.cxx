#include "Solver.hxx"
#include "SimulationProperties.hxx"
#include "NodePool.hxx"
#include "MultiDomain.hxx"
#include "FEMMath.hxx"

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

    namespace
    {
        std::tuple<SingleSolution, double, std::vector<double>>
          executeThermalSimulation(HygroThermFEM::MultiDomain & domain,
                                   const std::vector<double> & currentTemperature,
                                   const std::vector<double> & previousTimestepTemperature,
                                   const double dTime,
                                   size_t timestepIndex)
        {
            if(!domain.simulateThermal)
            {
                return std::make_tuple(
                  SingleSolution{currentTemperature, dTime}, 0.0, currentTemperature);
            }
            auto newTemperatureSolution =
              domain.thermalDomain.transient(previousTimestepTemperature, dTime, timestepIndex);
            auto newTemperatureError =
              HygroThermFEM::errorNorm(newTemperatureSolution.solution, currentTemperature);
            auto newCurrentTemperature = newTemperatureSolution.solution;
            return std::make_tuple(
              newTemperatureSolution, newTemperatureError, newCurrentTemperature);
        }

        std::tuple<SingleSolution, double, std::vector<double>>
          executeMoistureSimulation(HygroThermFEM::MultiDomain & domain,
                                    const std::vector<double> & currentHumidity,
                                    const std::vector<double> & previousTimestepHumidity,
                                    const double dTime,
                                    size_t timestepIndex)
        {
            if(!domain.simulateMoisture)
            {
                return std::make_tuple(
                  SingleSolution{currentHumidity, dTime}, 0.0, currentHumidity);
            }
            auto newHumiditySolution =
              domain.moistureDomain.transient(previousTimestepHumidity, dTime, timestepIndex);
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
                                       SingleSolution & temperatureSolution,
                                       SingleSolution & humiditySolution,
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
        SingleSolution temperatureSolution{previousTimestepTemperature, t_DTime};
        SingleSolution humiditySolution{previousTimestepHumidity, t_DTime};

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
                NodePool::Instance().updateNodeValues(
                  temperatureSolution.solution, BaseVariable::temperature, false);
                std::tie(humiditySolution, humidityError, currentHumidity) =
                  executeMoistureSimulation(
                    domain, currentHumidity, previousTimestepHumidity, dTime, timestepIndex);
            }

            if(domain.simulateThermal)
            {
                NodePool::Instance().updateNodeValues(
                  humiditySolution.solution, BaseVariable::humidity, false);
                std::tie(temperatureSolution, temperatureError, currentTemperature) =
                  executeThermalSimulation(
                    domain, currentTemperature, previousTimestepTemperature, dTime, timestepIndex);
            }

            ++currentIteration;
        }

        while((temperatureError > ConvergenceError && humidityError > ConvergenceError)
              || currentIteration > MaxIterations);

        NodePool::Instance().updateNodeValues(
          temperatureSolution.solution, BaseVariable::temperature, true);
        NodePool::Instance().updateNodeValues(
          humiditySolution.solution, BaseVariable::humidity, true);

        NodePool::Instance().updateNodeValues(currentHumidity, BaseVariable::humidity, true);
        NodePool::Instance().updateNodeValues(currentTemperature, BaseVariable::temperature, true);

        const auto waterContent = NodePool::Instance().properties(Variable::water);
        const auto liquidContent = NodePool::Instance().properties(Variable::liquid);
        const auto vaporContent = NodePool::Instance().properties(Variable::vapor);
        const auto iceContent = NodePool::Instance().properties(Variable::ice);

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
                humidity = domain.moistureDomain.steadyState();
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
                temperature = domain.thermalDomain.steadyState();
                temperatureError =
                  HygroThermFEM::errorNorm(temperature, previousTemperature);
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

        NodePool::Instance().updateNodeValues(humidity, BaseVariable::humidity, true);
        NodePool::Instance().updateNodeValues(temperature, BaseVariable::temperature, true);

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
}   // namespace HygroThermFEM