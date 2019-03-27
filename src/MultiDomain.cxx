#include <cmath>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "NodePool.hxx"

namespace HygroThermFEM
{
    // passing false to subdomains means that previous timestep values will not be automatically
    // updated. This mean that multidomain must update its values once solution converged.
    MultiDomain::MultiDomain(const bool performThermal, const bool performMoisture) :
        m_PerformThermal(performThermal),
        m_PerformMoisture(performMoisture)
    {}

    Solution MultiDomain::transient(std::vector<double> & temperature,
                                    std::vector<double> & humidity,
                                    const double t_DTime)
    {
        auto temperatureError = std::numeric_limits<double>::max();
        auto humidityError = std::numeric_limits<double>::max();
        auto currentTemperature = temperature;
        auto currentHumidity = humidity;
        double dTime = t_DTime;

        while(temperatureError > ConvergenceError || humidityError > ConvergenceError)
        {
            double dTimeThermal = t_DTime;
            double dTimeMoisture = t_DTime;
            SingleSolution temperatureSolution{temperature, t_DTime};
            SingleSolution humiditySolution{humidity, t_DTime};

            do   // Loop that performs adaptive timestep in case of convergence failure.
            {
                // do loop iterations need to make sure that both results are calculated for
                // identical timestep. This is part of adaptive timestep that program tries to
                // achieve in case when fail to converge.
                if(dTimeMoisture < dTimeThermal)
                {
                    dTimeThermal = dTimeMoisture;
                }

                if(dTimeMoisture > dTimeThermal)
                {
                    dTimeMoisture = dTimeThermal;
                }

                if(m_PerformThermal)
                {
                    m_ThermalDomain.updateNodeValues(
                      humiditySolution.solution, BaseVariable::humidity, false);
                    temperatureSolution = m_ThermalDomain.transient(temperature, dTimeThermal);
                    temperatureError = normError(temperatureSolution.solution, currentTemperature);
                    dTimeThermal = temperatureSolution.dTime;
                }
                else
                {
                    temperatureError = 0.0;
                }

                if(m_PerformMoisture)
                {
                    m_MoistureDomain.updateNodeValues(
                      temperatureSolution.solution, BaseVariable::temperature, false);
                    humiditySolution = m_MoistureDomain.transient(humidity, dTimeMoisture);
                    humidityError = normError(humiditySolution.solution, currentHumidity);
                    dTimeMoisture = humiditySolution.dTime;
                    dTime = dTimeThermal;
                }
                else
                {
                    humidityError = 0.0;
                }

            } while(dTimeThermal != dTimeMoisture);

            if(m_PerformThermal)
            {
                currentTemperature = temperatureSolution.solution;
            }

            if(m_PerformMoisture)
            {
                currentHumidity = humiditySolution.solution;
            }
        }

        m_ThermalDomain.updateNodeValues(currentTemperature, BaseVariable::temperature, true);
        m_ThermalDomain.updateNodeValues(currentHumidity, BaseVariable::humidity, true);
        m_MoistureDomain.updateNodeValues(currentTemperature, BaseVariable::temperature, true);
        m_MoistureDomain.updateNodeValues(currentHumidity, BaseVariable::humidity, true);

        NodePool::Instance().updateNodeValues(currentHumidity, BaseVariable::humidity);
        NodePool::Instance().updateNodeValues(currentTemperature, BaseVariable::temperature);

        const auto waterContent = NodePool::Instance().properties(Variable::water);
        const auto liquidContent = NodePool::Instance().properties(Variable::liquid);
        const auto vaporContent = NodePool::Instance().properties(Variable::vapor);
        const auto iceContent = NodePool::Instance().properties(Variable::ice);

        return Solution{dTime,
                        currentTemperature,
                        currentHumidity,
                        waterContent,
                        liquidContent,
                        vaporContent,
                        iceContent};
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
                                              const double t_AirTemperature,
                                              const double t_ConvectionCoefficient,
                                              const double t_Humidity)
    {
        m_ThermalDomain.createConvectionBCFixedHc(
          index1, index2, t_AirTemperature, t_ConvectionCoefficient, t_Humidity, m_PerformMoisture);

        m_MoistureDomain.createMoistureBCFixedHc(
          index1, index2, t_AirTemperature, t_ConvectionCoefficient, t_Humidity);
    }

    void MultiDomain::createMoistureBCVariableHc(const size_t index1,
                                                 const size_t index2,
                                                 const double t_AirTemperature,
                                                 const double t_Humidity)
    {
        m_ThermalDomain.createConvectionBCVariableHc(index1, index2, t_AirTemperature, t_Humidity, m_PerformMoisture);

        m_MoistureDomain.createMoistureBCVariableHc(index1, index2, t_Humidity, t_AirTemperature);
    }

    void MultiDomain::createTemperatureBC(const size_t index1,
                                          const size_t index2,
                                          const double t_Temp1,
                                          const double t_Temp2)
    {
        m_ThermalDomain.createTemperatureBC(index1, index2, t_Temp1, t_Temp2);
    }

    void MultiDomain::createTemperatureBC(const size_t index1,
                                          const size_t index2,
                                          const double t_Temp)
    {
        m_ThermalDomain.createTemperatureBC(index1, index2, t_Temp);
    }

    void MultiDomain::createBlackBodyRadiationBC(const size_t index1,
                                                 const size_t index2,
                                                 const double t_Emissivity,
                                                 const double t_RadiationTemperature)
    {
        m_ThermalDomain.createBlackBodyRadiationBC(
          index1, index2, t_Emissivity, t_RadiationTemperature);
    }

    void MultiDomain::createSimplifiedRadiationBC(const size_t index1,
                                                  const size_t index2,
                                                  const double t_RadiationCoefficient,
                                                  const double t_RadiationTemperature)
    {
        m_ThermalDomain.createSimplifiedRadiationBC(
          index1, index2, t_RadiationCoefficient, t_RadiationTemperature);
    }

    double MultiDomain::normError(const std::vector<double> & vec1,
                                  const std::vector<double> & vec2)
    {
        const auto norm1 = norm(vec1);
        const auto norm2 = norm(vec2);

        return std::abs(norm1 - norm2) / norm1;
    }

    std::vector<double> MultiDomain::property(Variable property)
    {
        return NodePool::Instance().properties(property);
    }

    Solution::Solution(const double dtime,
                       const std::vector<double> & temperature,
                       const std::vector<double> & humidity,
                       const std::vector<double> & waterContent,
                       const std::vector<double> & liquidWaterContent,
                       const std::vector<double> & vaporContent,
                       const std::vector<double> & iceContent) :
        dTime(dtime),
        temperature(temperature),
        humidity(humidity),
        waterContent(waterContent),
        liquidWaterContent(liquidWaterContent),
        vaporContent(vaporContent),
        iceContent(iceContent)
    {}
}   // namespace HygroThermFEM
