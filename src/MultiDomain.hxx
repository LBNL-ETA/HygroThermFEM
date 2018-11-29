#pragma once

#include "Domains.hxx"

namespace MoisThermFEM
{
    struct Solution
    {
        Solution(const std::vector<double> & temperature,
                 const std::vector<double> & humidity,
                 const std::vector<double> & waterContent,
                 const std::vector<double> & liquidWaterContent,
                 const std::vector<double> & vaporContent,
                 const std::vector<double> & iceContent);

        std::vector<double> temperature;
        std::vector<double> humidity;
        std::vector<double> waterContent;
        std::vector<double> liquidWaterContent;
        std::vector<double> vaporContent;
        std::vector<double> iceContent;
    };

    class MultiDomain
    {
    public:
        MultiDomain();

        /// Calculates next timestep value from current values
        Solution transient(std::vector<double> & temperature,
                           std::vector<double> & humidity,
                           const double t_DTime);

        std::vector<double> property(Variable property) const;

        void createElement(const size_t index1,
                           const size_t index2,
                           const size_t index3,
                           const size_t index4,
                           const std::string & materialName);

        void createConvectionBC(const size_t index1,
                                const size_t index2,
                                const double t_ConvectionCoefficient,
                                const double t_AirTemperature,
                                const double t_Humidity);

        void createTemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_Temp1,
                                 const double t_Temp2);

        void createTemperatureBC(const size_t index1, const size_t index2, const double t_Temp);

        void createBlackBodyRadiationBC(const size_t index1,
                                        const size_t index2,
                                        const double t_Emissivity,
                                        const double t_RadiationTemperature);

    private:
        static double normError(const std::vector<double> & vec1, const std::vector<double> & vec2);

        ThermalDomain m_ThermalDomain;
        MoistureDomain m_MoistureDomain;
    };

}   // namespace MoisThermFEM
