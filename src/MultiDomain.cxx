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
    // passing false to subdomains means that previous timestep values will not be automatically
    // updated. This mean that multidomain must update its values once solution converged.
    MultiDomain::MultiDomain(const bool performThermal, const bool performMoisture) :
        m_SimulateThermal(performThermal), m_SimulateMoisture(performMoisture)
    {}

    Solution MultiDomain::transient(const std::vector<double> & temperature,
                                    const std::vector<double> & humidity,
                                    const double t_DTime,
                                    size_t timestepIndex)
    {
        const auto ConvergenceError{SimulationProperties::Instance().errorTolerance()};
        auto temperatureError{std::numeric_limits<double>::max()};
        auto humidityError{std::numeric_limits<double>::max()};
        auto currentTemperature{temperature};
        auto currentHumidity{humidity};
        double dTime{t_DTime};
        SingleSolution temperatureSolution{temperature, t_DTime};
        SingleSolution humiditySolution{humidity, t_DTime};

        const auto MaxIterations = SimulationProperties::Instance().maxNumberOfIterations();

        size_t currentIteration{0};

        // Note that temperature and humidity are solved separately first and then updated with new
        // data for next iteration.
        do
        {
            executeTransientIteration(temperature,
                                      humidity,
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

            if(m_SimulateMoisture)
            {
                NodePool::Instance().updateNodeValues(
                  temperatureSolution.solution, BaseVariable::temperature, false);
                std::tie(humiditySolution, humidityError, currentHumidity) =
                  executeMoistureSimulation(currentHumidity, humidity, dTime, timestepIndex);
            }

            if(m_SimulateThermal)
            {
                NodePool::Instance().updateNodeValues(
                  humiditySolution.solution, BaseVariable::humidity, false);
                std::tie(temperatureSolution, temperatureError, currentTemperature) =
                  executeThermalSimulation(currentTemperature, temperature, dTime, timestepIndex);
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

        const auto heatFlux = m_ThermalDomain.flux();
        const auto waterFlux = m_MoistureDomain.flux();

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
            if(m_SimulateMoisture)
            {
                humidity = m_MoistureDomain.steadyState();
                humidityError = normError(humidity, previousHumidity);
                previousHumidity = humidity;
                NodePool::Instance().updateNodeValues(humidity, BaseVariable::humidity);
            }
            else
            {
                humidityError = 0;
            }
            if(m_SimulateThermal)
            {
                temperature = m_ThermalDomain.steadyState();
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

        const auto heatFlux = m_ThermalDomain.flux();
        const auto waterFlux = m_MoistureDomain.flux();

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

    void MultiDomain::performMoistureSimulation(const bool val)
    {
        m_SimulateMoisture = val;
    }

    void MultiDomain::performThermalSimulation(const bool val)
    {
        m_SimulateThermal = val;
    }

    void MultiDomain::createElement(const size_t index1,
                                    const size_t index2,
                                    const size_t index3,
                                    const size_t index4,
                                    const std::string & materialName)
    {
        m_ThermalDomain.createElement(index1, index2, index3, index4, materialName);
        m_MoistureDomain.createElement(index1, index2, index3, index4, materialName);
    }

    void MultiDomain::createBC_FixedHc(const size_t index1,
                                       const size_t index2,
                                       const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        m_ThermalDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients, m_SimulateMoisture);

        m_MoistureDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients);
    }

    void MultiDomain::createBC_FixedHc(
      size_t index1,
      size_t index2,
      const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        m_ThermalDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients, m_SimulateMoisture);
        m_MoistureDomain.createBC_FixedHc(index1, index2, fixedBchcCoefficients);
    }

    void MultiDomain::createBC_TARPHc(size_t index1,
                                      size_t index2,
                                      const TARPCoefficients & varHCCoeff,
                                      double surfaceTilt)
    {
        m_ThermalDomain.createBC_TARPHc(
          index1, index2, varHCCoeff, surfaceTilt, m_SimulateMoisture);

        m_MoistureDomain.createBC_TARPHc(index1, index2, varHCCoeff, surfaceTilt);
    }

    void MultiDomain::createBC_TARPHc(size_t index1,
                                      size_t index2,
                                      const std::vector<TARPCoefficients> & varHCCoeff,
                                      double surfaceTilt)
    {
        m_ThermalDomain.createBC_TARPHc(
          index1, index2, varHCCoeff, surfaceTilt, m_SimulateMoisture);

        m_MoistureDomain.createBC_TARPHc(index1, index2, varHCCoeff, surfaceTilt);
    }

    void MultiDomain::createBC_ASHRAEInsideHc(size_t index1,
                                              size_t index2,
                                              const ASHRAEInsideCoefficients & coeff,
                                              double surfaceHeight,
                                              double surfaceTilt)
    {
        m_ThermalDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt, m_SimulateMoisture);
        m_MoistureDomain.createBC_ASHRAEInsideHc(index1, index2, coeff, surfaceHeight, surfaceTilt);
    }

    void MultiDomain::createBC_ASHRAEInsideHc(size_t index1,
                                              size_t index2,
                                              const std::vector<ASHRAEInsideCoefficients> & coeff,
                                              double surfaceHeight,
                                              double surfaceTilt)
    {
        m_ThermalDomain.createBC_ASHRAEInsideHc(
          index1, index2, coeff, surfaceHeight, surfaceTilt, m_SimulateMoisture);
        m_MoistureDomain.createBC_ASHRAEInsideHc(index1, index2, coeff, surfaceHeight, surfaceTilt);
    }

    void MultiDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                               size_t index2,
                                               const ASHRAEOutsideCoefficients & coeff)
    {
        m_ThermalDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff, m_SimulateMoisture);
        m_MoistureDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff);
    }

    void MultiDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                               size_t index2,
                                               const std::vector<ASHRAEOutsideCoefficients> & coeff)
    {
        m_ThermalDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff, m_SimulateMoisture);
        m_MoistureDomain.createBC_ASHRAEOutsideHc(index1, index2, coeff);
    }

    void MultiDomain::createBC_YazdanianKlemsHc(size_t index1,
                                                size_t index2,
                                                const YazdanianKlemsCoefficients & coeff)
    {
        m_ThermalDomain.createBC_YazdanianKlemsHc(index1, index2, coeff, m_SimulateMoisture);
        m_MoistureDomain.createBC_YazdanianKlemsHc(index1, index2, coeff);
    }

    void MultiDomain::createBC_YazdanianKlemsHc(
      size_t index1, size_t index2, const std::vector<YazdanianKlemsCoefficients> & coeff)
    {
        m_ThermalDomain.createBC_YazdanianKlemsHc(index1, index2, coeff, m_SimulateMoisture);
        m_MoistureDomain.createBC_YazdanianKlemsHc(index1, index2, coeff);
    }

    void
      MultiDomain::createBC_KimuraHc(size_t index1, size_t index2, const KimuraCoefficients & coeff)
    {
        m_ThermalDomain.createBC_KimuraHc(index1, index2, coeff, m_SimulateMoisture);
        m_MoistureDomain.createBC_KimuraHc(index1, index2, coeff);
    }

    void MultiDomain::createBC_KimuraHc(size_t index1,
                                        size_t index2,
                                        const std::vector<KimuraCoefficients> & coeff)
    {
        m_ThermalDomain.createBC_KimuraHc(index1, index2, coeff, m_SimulateMoisture);
        m_MoistureDomain.createBC_KimuraHc(index1, index2, coeff);
    }

    void MultiDomain::createBC_FixedTemperature(const size_t index1,
                                                const size_t index2,
                                                const double t_Temp1,
                                                const double t_Temp2)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, t_Temp1, t_Temp2);
    }

    void MultiDomain::createBC_FixedTemperature(size_t index1,
                                                size_t index2,
                                                const std::vector<ConstantBCTemperatures> & temp)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, temp);
    }

    void MultiDomain::createBC_FixedTemperature(const size_t index1,
                                                const size_t index2,
                                                const double t_Temp)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, t_Temp);
    }

    void
      MultiDomain::createBC_FixedTemperature(size_t index1, size_t index2, std::vector<double> temp)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, std::move(temp));
    }

    void MultiDomain::createBC_FixedTemperatureAndHumidity(size_t index1,
                                                           size_t index2,
                                                           const TemperatureAndHumidity & values)
    {
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, values.Temperature);
        m_MoistureDomain.createBC_FixedHumidity(index1, index2, values);
    }

    void MultiDomain::createBC_FixedTemperatureAndHumidity(
      size_t index1, size_t index2, const std::vector<TemperatureAndHumidity> & values)
    {
        std::vector<double> temperatures(values.size());
        for(size_t i = 0u; i < values.size(); ++i)
        {
            temperatures[i] = values[i].Temperature;
        }
        m_ThermalDomain.createBC_FixedTemperature(index1, index2, temperatures);
        m_MoistureDomain.createBC_FixedHumidity(index1, index2, values);
    }

    void MultiDomain::createBC_FixedHeatFlux(size_t index1, size_t index2, double t_Flux)
    {
        m_ThermalDomain.createBC_FixedFlux(index1, index2, t_Flux);
    }

    void
      MultiDomain::createBC_FixedHeatFlux(size_t index1, size_t index2, std::vector<double> t_Flux)
    {
        m_ThermalDomain.createBC_FixedFlux(index1, index2, t_Flux);
    }

    void MultiDomain::createBC_BlackBodyRadiation(size_t index1,
                                                  size_t index2,
                                                  double t_Emissivity,
                                                  double t_RadiationTemperature)
    {
        m_ThermalDomain.createBC_BlackBodyRadiation(
          index1, index2, t_Emissivity, t_RadiationTemperature);
    }

    void MultiDomain::createBC_BlackBodyRadiation(
      size_t index1, size_t index2, const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        m_ThermalDomain.createBC_BlackBodyRadiation(index1, index2, radCoeffs);
    }

    void MultiDomain::createBC_LinearizedRadiation(
      const size_t index1,
      const size_t index2,
      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        m_ThermalDomain.createBC_LinearizedRadiation(index1, index2, linearRadBC);
    }

    void MultiDomain::createBC_LinearizedRadiation(
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC)
    {
        m_ThermalDomain.createBC_LinearizedRadiation(index1, index2, linearRadBC);
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
        m_ThermalDomain.setGravityVector(gravityVector);
    }

    void MultiDomain::subscribeThermal(Timesteps::TimestepObserver * observer)
    {
        m_ThermalDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeThermal(Timesteps::TimestepObserver * observer)
    {
        m_ThermalDomain.unsubscribe(observer);
    }

    void MultiDomain::subscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        m_MoistureDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        m_MoistureDomain.unsubscribe(observer);
    }

    bool MultiDomain::isMoistureSimulationON() const
    {
        return m_SimulateMoisture;
    }

    bool MultiDomain::isThermalSimulationON() const
    {
        return m_SimulateThermal;
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
        m_ThermalDomain.clearModel();
        m_MoistureDomain.clearModel();
    }
    std::tuple<SingleSolution, double, std::vector<double>>
      MultiDomain::executeThermalSimulation(const std::vector<double> & currentTemperature,
                                            const std::vector<double> & temperature,
                                            const double dTime,
                                            size_t timestepIndex)
    {
        if(!m_SimulateThermal)
        {
            return std::make_tuple(SingleSolution{currentTemperature, dTime}, 0.0, currentTemperature);
        }
        auto newTemperatureSolution = m_ThermalDomain.transient(temperature, dTime, timestepIndex);
        auto newTemperatureError = normError(newTemperatureSolution.solution, currentTemperature);
        auto newCurrentTemperature = newTemperatureSolution.solution;
        return std::make_tuple(newTemperatureSolution, newTemperatureError, newCurrentTemperature);
    }

    std::tuple<SingleSolution, double, std::vector<double>>
      MultiDomain::executeMoistureSimulation(const std::vector<double> & currentHumidity,
                                             const std::vector<double> & humidity,
                                             const double dTime,
                                             size_t timestepIndex)
    {
        if(!m_SimulateMoisture)
        {
            return std::make_tuple(SingleSolution{currentHumidity, dTime}, 0.0, currentHumidity);
        }
        auto newHumiditySolution = m_MoistureDomain.transient(humidity, dTime, timestepIndex);
        auto newHumidityError = normError(newHumiditySolution.solution, currentHumidity);
        auto newCurrentHumidity = newHumiditySolution.solution;
        return std::make_tuple(newHumiditySolution, newHumidityError, newCurrentHumidity);
    }

    void MultiDomain::executeTransientIteration(const std::vector<double> & temperature,
                                                const std::vector<double> & humidity,
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
              executeMoistureSimulation(currentHumidity, humidity, dTime, timestepIndex);
            std::tie(temperatureSolution, temperatureError, currentTemperature) =
              executeThermalSimulation(currentTemperature, temperature, dTime, timestepIndex);
            ++localIterCounter;
        }
    }

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
