#include <cmath>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "NodePool.hxx"

namespace MoisThermFEM
{
    MultiDomain::MultiDomain() : m_ThermalDomain(), m_MoistureDomain()
    {}

    Solution MultiDomain::transient(std::vector<double> & temperature,
                                    std::vector<double> & humidity,
                                    const double t_DTime)
    {
        auto temperatureError = std::numeric_limits<double>::max();
        auto humidityError = std::numeric_limits<double>::max();
        auto currentTemperature = temperature;
        auto currentHumidity = humidity;

        while(temperatureError > ConvergenceError || humidityError > ConvergenceError)
        {
            m_ThermalDomain.updateNodeValues(currentHumidity, BaseVariable::humidity);
            const auto temperatureSolution = m_ThermalDomain.transient(temperature, t_DTime);
            temperatureError = normError(temperatureSolution.solution, currentTemperature);

            m_MoistureDomain.updateNodeValues(currentTemperature, BaseVariable::temperature);
            const auto humiditySolution = m_MoistureDomain.transient(humidity, t_DTime);
            humidityError = normError(humiditySolution.solution, currentHumidity);

            currentHumidity = humiditySolution.solution;
            currentTemperature = temperatureSolution.solution;
        }

        NodePool::Instance().updateNodeValues(currentHumidity, BaseVariable::humidity);
        NodePool::Instance().updateNodeValues(currentTemperature, BaseVariable::temperature);

        const auto waterContent = NodePool::Instance().properties(Variable::water);
        const auto liquidContent = NodePool::Instance().properties(Variable::liquid);
        const auto vaporContent = NodePool::Instance().properties(Variable::vapor);
        const auto iceContent = NodePool::Instance().properties(Variable::ice);

        return Solution{currentTemperature,
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

    void MultiDomain::createConvectionBC(const size_t index1,
                                         const size_t index2,
                                         const double t_ConvectionCoefficient,
                                         const double t_AirTemperature,
                                         const double t_Humidity)
    {
        m_ThermalDomain.createConvectionBC(
          index1, index2, t_ConvectionCoefficient, t_AirTemperature);

		m_MoistureDomain.createMoistureBC(
			index1, index2, t_Humidity, t_AirTemperature );
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

    double MultiDomain::normError(const std::vector<double> & vec1,
                                  const std::vector<double> & vec2)
    {
        const auto norm1 = norm(vec1);
        const auto norm2 = norm(vec2);

        return std::abs(norm1 - norm2) / norm1;
    }

    std::vector<double> MultiDomain::property(Variable property) const
    {
        return NodePool::Instance().properties(property);
    }

    Solution::Solution(const std::vector<double> & temperature,
                       const std::vector<double> & humidity,
                       const std::vector<double> & waterContent,
                       const std::vector<double> & liquidWaterContent,
                       const std::vector<double> & vaporContent,
                       const std::vector<double> & iceContent) :
        temperature(temperature),
        humidity(humidity),
        waterContent(waterContent),
        liquidWaterContent(liquidWaterContent),
        vaporContent(vaporContent),
        iceContent(iceContent)
    {}
}   // namespace MoisThermFEM
