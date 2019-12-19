#include <cmath>
#include <limits>
#include <utility>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "NodePool.hxx"
#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    // passing false to subdomains means that previous timestep values will not be automatically
    // updated. This mean that multidomain must update its values once solution converged.
    MultiDomain::MultiDomain(const bool performThermal, const bool performMoisture) :
        m_PerformThermal(performThermal), m_PerformMoisture(performMoisture)
    {}

    Solution MultiDomain::transient(std::vector<double> & temperature,
                                    std::vector<double> & humidity,
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
            size_t localIterCounter{0};
            while(humidityError > ConvergenceError && localIterCounter <= MaxIterations
                  && temperatureError > ConvergenceError && localIterCounter <= MaxIterations)
            {
                if(m_PerformMoisture)
                {
                    humiditySolution = m_MoistureDomain.transient(humidity, dTime, timestepIndex);
                    humidityError = normError(humiditySolution.solution, currentHumidity);
                    currentHumidity = humiditySolution.solution;
                    ++localIterCounter;
                }
                else
                {
                    humidityError = 0;
                }
                if(m_PerformThermal)
                {
                    temperatureSolution =
                            m_ThermalDomain.transient(temperature, dTime, timestepIndex);
                    temperatureError = normError(temperatureSolution.solution, currentTemperature);
                    currentTemperature = temperatureSolution.solution;
                    ++localIterCounter;
                }
                else
                {
                    temperatureError = 0;
                }
            }


            if(m_PerformMoisture)
            {
                NodePool::Instance().updateNodeValues(
                        temperatureSolution.solution, BaseVariable::temperature, false);
                humiditySolution = m_MoistureDomain.transient(humidity, dTime, timestepIndex);
                humidityError = normError(humiditySolution.solution, currentHumidity);
                currentHumidity = humiditySolution.solution;
            }

            if(m_PerformThermal)
            {
                NodePool::Instance().updateNodeValues(
                        humiditySolution.solution, BaseVariable::humidity, false);
                temperatureSolution = m_ThermalDomain.transient(temperature, dTime, timestepIndex);
                temperatureError = normError(temperatureSolution.solution, currentTemperature);
                currentTemperature = temperatureSolution.solution;
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
            if(m_PerformMoisture)
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
            if(m_PerformThermal)
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
        m_PerformMoisture = val;
    }

    void MultiDomain::performThermalSimulation(const bool val)
    {
        m_PerformThermal = val;
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

    void MultiDomain::createMoistureBCFixedHc(const size_t index1,
                                              const size_t index2,
                                              const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        m_ThermalDomain.createConvectionBCFixedHc(
          index1, index2, fixedBchcCoefficients, m_PerformMoisture);

        m_MoistureDomain.createMoistureBCFixedHc(index1, index2, fixedBchcCoefficients);
    }

    void MultiDomain::createMoistureBCFixedHc(
      size_t index1,
      size_t index2,
      const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        m_ThermalDomain.createConvectionBCFixedHc(
          index1, index2, fixedBchcCoefficients, m_PerformMoisture);
        m_MoistureDomain.createMoistureBCFixedHc(index1, index2, fixedBchcCoefficients);
    }

    void MultiDomain::createMoistureBCVariableHc(const size_t index1,
                                                 const size_t index2,
                                                 const VariableBCHCCoefficients & varHCCoeff)
    {
        m_ThermalDomain.createConvectionBCVariableHc(index1, index2, varHCCoeff, m_PerformMoisture);

        m_MoistureDomain.createMoistureBCVariableHc(index1, index2, varHCCoeff);
    }

    void MultiDomain::createMoistureBCVariableHc(
      size_t index1, size_t index2, const std::vector<VariableBCHCCoefficients> & varHCCoeff)
    {
        m_ThermalDomain.createConvectionBCVariableHc(index1, index2, varHCCoeff, m_PerformMoisture);

        m_MoistureDomain.createMoistureBCVariableHc(index1, index2, varHCCoeff);
    }

    void MultiDomain::createTemperatureBC(const size_t index1,
                                          const size_t index2,
                                          const double t_Temp1,
                                          const double t_Temp2)
    {
        m_ThermalDomain.createTemperatureBC(index1, index2, t_Temp1, t_Temp2);
    }

    void MultiDomain::createTemperatureBC(size_t index1,
                                          size_t index2,
                                          const std::vector<ConstantBCTemperatures> & temp)
    {
        m_ThermalDomain.createTemperatureBC(index1, index2, temp);
    }

    void MultiDomain::createTemperatureBC(const size_t index1,
                                          const size_t index2,
                                          const double t_Temp)
    {
        m_ThermalDomain.createTemperatureBC(index1, index2, t_Temp);
    }

    void MultiDomain::createTemperatureBC(size_t index1, size_t index2, std::vector<double> temp)
    {
        m_ThermalDomain.createTemperatureBC(index1, index2, std::move(temp));
    }

    void MultiDomain::createBlackBodyRadiationBC(const size_t index1,
                                                 const size_t index2,
                                                 const double t_Emissivity,
                                                 const double t_RadiationTemperature)
    {
        m_ThermalDomain.createBlackBodyRadiationBC(
          index1, index2, t_Emissivity, t_RadiationTemperature);
    }

    void MultiDomain::createBlackBodyRadiationBC(
      size_t index1, size_t index2, const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        m_ThermalDomain.createBlackBodyRadiationBC(index1, index2, radCoeffs);
    }

    void MultiDomain::createLinearizedRadiationBC(
      const size_t index1,
      const size_t index2,
      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        m_ThermalDomain.createLinearizedRadiationBC(index1, index2, linearRadBC);
    }

    void MultiDomain::createLinearizedRadiationBC(
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC)
    {
        m_ThermalDomain.createLinearizedRadiationBC(index1, index2, linearRadBC);
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
