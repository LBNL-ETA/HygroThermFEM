#include <cmath>
#include <limits>
#include <utility>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "NodePool.hxx"
#include "SimulationProperties.hxx"
#include "MaterialDataChecker.hxx"

namespace HygroThermFEM
{
    Solution MultiDomain::transient(const std::vector<double> & previousTimestepTemperature,
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
            executeTransientIteration(previousTimestepTemperature,
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

            if(simulateMoisture)
            {
                NodePool::Instance().updateNodeValues(
                  temperatureSolution.solution, BaseVariable::temperature, false);
                std::tie(humiditySolution, humidityError, currentHumidity) =
                  executeMoistureSimulation(
                    currentHumidity, previousTimestepHumidity, dTime, timestepIndex);
            }

            if(simulateThermal)
            {
                NodePool::Instance().updateNodeValues(
                  humiditySolution.solution, BaseVariable::humidity, false);
                std::tie(temperatureSolution, temperatureError, currentTemperature) =
                  executeThermalSimulation(
                    currentTemperature, previousTimestepTemperature, dTime, timestepIndex);
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

        const auto heatFlux = thermalDomain.flux();
        const auto waterFlux = moistureDomain.flux();

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

    Solution MultiDomain::steadyState()
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
            if(simulateMoisture)
            {
                humidity = moistureDomain.steadyState();
                humidityError = normError(humidity, previousHumidity);
                previousHumidity = humidity;
                NodePool::Instance().updateNodeValues(humidity, BaseVariable::humidity);
            }
            else
            {
                humidityError = 0;
            }
            if(simulateThermal)
            {
                temperature = thermalDomain.steadyState();
                temperatureError = normError(temperature, previousTemperature);
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

        const auto waterContent = NodePool::Instance().properties(Variable::water);
        const auto liquidContent = NodePool::Instance().properties(Variable::liquid);
        const auto vaporContent = NodePool::Instance().properties(Variable::vapor);
        const auto iceContent = NodePool::Instance().properties(Variable::ice);

        const auto heatFlux = thermalDomain.flux();
        const auto waterFlux = moistureDomain.flux();

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

    void MultiDomain::createElement(const size_t index1,
                                    const size_t index2,
                                    const size_t index3,
                                    const size_t index4,
                                    const std::string & materialName)
    {
        thermalDomain.createElement(index1, index2, index3, index4, materialName);
        moistureDomain.createElement(index1, index2, index3, index4, materialName);
    }

    double MultiDomain::normError(const std::vector<double> & vec1,
                                  const std::vector<double> & vec2)
    {
        auto norm1 = norm(vec1);
        auto norm2 = norm(vec2);
        if(norm1 == 0)
        {
            norm1 = 1e-10;
            if(norm2 == 0)
            {
                norm2 = norm1;
            }
        }

        return std::abs(norm1 - norm2) / norm1;
    }

    std::vector<double> MultiDomain::property(Variable property)
    {
        return NodePool::Instance().properties(property);
    }

    void MultiDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        thermalDomain.setGravityVector(gravityVector);
    }

    void MultiDomain::subscribeThermal(Timesteps::TimestepObserver * observer)
    {
        thermalDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeThermal(Timesteps::TimestepObserver * observer)
    {
        thermalDomain.unsubscribe(observer);
    }

    void MultiDomain::subscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        moistureDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        moistureDomain.unsubscribe(observer);
    }

    MaterialsErrorCheckVector MultiDomain::checkMaterialsForTransientSimulation() const
    {
        MaterialDataChecker dataChecker{*this};
        return dataChecker.checkMaterialProperties(true);
    }

    MaterialsErrorCheckVector MultiDomain::checkMaterialsForSteadyStateSimulation() const
    {
        MaterialDataChecker dataChecker{*this};
        return dataChecker.checkMaterialProperties(false);
    }

    MaterialsErrorCheckVector
      MultiDomain::checkForMaterialsValidity(const SimulationType simulationType) const
    {
        MaterialsErrorCheckVector result;
        switch(simulationType)
        {
            case SimulationType::SteadyState:
                result = checkMaterialsForSteadyStateSimulation();
                break;
            case SimulationType::Transient:
                result = checkMaterialsForTransientSimulation();
                break;
            default:
                throw std::runtime_error("Incorrect selection of simulation type.");
        }
        return result;
    }

    void MultiDomain::clearModel()
    {
        thermalDomain.clearModel();
        moistureDomain.clearModel();
    }
    std::tuple<SingleSolution, double, std::vector<double>>
      MultiDomain::executeThermalSimulation(const std::vector<double> & currentTemperature,
                                            const std::vector<double> & previousTimestepTemperature,
                                            const double dTime,
                                            size_t timestepIndex)
    {
        if(!simulateThermal)
        {
            return std::make_tuple(
              SingleSolution{currentTemperature, dTime}, 0.0, currentTemperature);
        }
        auto newTemperatureSolution =
          thermalDomain.transient(previousTimestepTemperature, dTime, timestepIndex);
        auto newTemperatureError = normError(newTemperatureSolution.solution, currentTemperature);
        auto newCurrentTemperature = newTemperatureSolution.solution;
        return std::make_tuple(newTemperatureSolution, newTemperatureError, newCurrentTemperature);
    }

    std::tuple<SingleSolution, double, std::vector<double>>
      MultiDomain::executeMoistureSimulation(const std::vector<double> & currentHumidity,
                                             const std::vector<double> & previousTimestepHumidity,
                                             const double dTime,
                                             size_t timestepIndex)
    {
        if(!simulateMoisture)
        {
            return std::make_tuple(SingleSolution{currentHumidity, dTime}, 0.0, currentHumidity);
        }
        auto newHumiditySolution =
          moistureDomain.transient(previousTimestepHumidity, dTime, timestepIndex);
        auto newHumidityError = normError(newHumiditySolution.solution, currentHumidity);
        auto newCurrentHumidity = newHumiditySolution.solution;
        return std::make_tuple(newHumiditySolution, newHumidityError, newCurrentHumidity);
    }

    void MultiDomain::executeTransientIteration(
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
            std::tie(humiditySolution, humidityError, currentHumidity) = executeMoistureSimulation(
              currentHumidity, previousTimestepHumidity, dTime, timestepIndex);
            std::tie(temperatureSolution, temperatureError, currentTemperature) =
              executeThermalSimulation(
                currentTemperature, previousTimestepTemperature, dTime, timestepIndex);
            ++localIterCounter;
        }
    }

    MultiDomain::MultiDomain(bool performThermal, bool performMoisture) :
        simulateThermal(performThermal), simulateMoisture(performMoisture)
    {}

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
}   // namespace HygroThermFEM
