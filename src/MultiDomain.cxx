#include <cmath>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "NodePool.hxx"

namespace MoisThermFEM
{
    MultiDomain::MultiDomain() :
        m_ThermalDomain(),
        m_MoistureDomain()
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
            m_ThermalDomain.updateNodeValues(currentHumidity, StateProperty::humidity);
            const auto temperatureSolution = m_ThermalDomain.transient(temperature, t_DTime);
            temperatureError = normError(temperatureSolution, currentTemperature);

            m_MoistureDomain.updateNodeValues(currentTemperature, StateProperty::temperature);
            const auto humiditySolution = m_MoistureDomain.transient(humidity, t_DTime);
            humidityError = normError(humiditySolution, currentHumidity);

            currentHumidity = humiditySolution;
            currentTemperature = temperatureSolution;
        }

		NodePool::Instance().updateNodeValues(currentHumidity, StateProperty::humidity);
		NodePool::Instance().updateNodeValues(currentTemperature, StateProperty::temperature);

        const auto waterContent = NodePool::Instance().properties(Property::water);
        const auto liquidContent = NodePool::Instance().properties(Property::liquid);
        const auto vaporContent = NodePool::Instance().properties(Property::vapor);
        const auto iceContent = NodePool::Instance().properties(Property::ice);

        return Solution{currentTemperature, currentHumidity, waterContent,
						liquidContent, vaporContent, iceContent};
    }

    void MultiDomain::createElement(Node2D &t_Node1,
                                    Node2D &t_Node2,
                                    Node2D &t_Node3,
                                    Node2D &t_Node4,
                                    const Material &mat)
    {
        m_ThermalDomain.createElement(t_Node1, t_Node2, t_Node3, t_Node4, mat);
        m_MoistureDomain.createElement(t_Node1, t_Node2, t_Node3, t_Node4, mat);
    }

    void MultiDomain::createConvectionBC(Node2D &t_Node1,
                                         Node2D &t_Node2,
                                         const double t_ConvectionCoefficient,
                                         const double t_AirTemperature,
                                         const double t_Humidity)
    {
        m_ThermalDomain.createConvectionBC(
          t_Node1, t_Node2, t_ConvectionCoefficient, t_AirTemperature);

        m_MoistureDomain.createMoistureBC(
          t_Node1, t_Node2, t_ConvectionCoefficient, t_Humidity, t_AirTemperature);
    }

    void MultiDomain::createTemperatureBC(Node2D & t_Node1,
                                          Node2D & t_Node2,
                                          const double t_Temp1,
                                          const double t_Temp2)
    {
        m_ThermalDomain.createTemperatureBC(t_Node1, t_Node2, t_Temp1, t_Temp2);
    }

    void MultiDomain::createTemperatureBC(Node2D & t_Node1, Node2D & t_Node2, const double t_Temp)
    {
        m_ThermalDomain.createTemperatureBC(t_Node1, t_Node2, t_Temp);
    }

    void MultiDomain::createBlackBodyRadiationBC(Node2D &t_Node1,
                                                 Node2D &t_Node2,
                                                 const double t_Emissivity,
                                                 const double t_RadiationTemperature)
    {
        m_ThermalDomain.createBlackBodyRadiationBC(
          t_Node1, t_Node2, t_Emissivity, t_RadiationTemperature);
    }

    double MultiDomain::normError(const std::vector<double> & vec1,
                                  const std::vector<double> & vec2)
    {
        const auto norm1 = norm(vec1);
        const auto norm2 = norm(vec2);

        return std::abs(norm1 - norm2) / norm1;
    }

    std::vector<double> MultiDomain::property(Property property) const
    {
        return NodePool::Instance().properties(property);
    }

    Solution::Solution(const std::vector<double> & temperature,
                       const std::vector<double> & humidity,
                       const std::vector<double> & waterContent,
					   const std::vector<double> & liquidWaterContent,
					   const std::vector<double> & vaporContent,
					   const std::vector<double> & iceContent
                       ) :
        temperature(temperature),
        humidity(humidity),
        waterContent(waterContent),
        liquidWaterContent(liquidWaterContent),
        vaporContent(vaporContent),
        iceContent(iceContent)
    {}
}   // namespace MoisThermFEM
