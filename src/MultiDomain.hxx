#pragma once

#include "Domains.hxx"

namespace HygroThermFEM
{
    struct Solution
    {
        Solution(double dtime,
                 const std::vector<double> & temperature,
                 const std::vector<double> & humidity,
                 const std::vector<double> & waterContent,
                 const std::vector<double> & liquidWaterContent,
                 const std::vector<double> & vaporContent,
                 const std::vector<double> & iceContent);

        double dTime;
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
                           double t_DTime);

        static std::vector<double> property(Variable property);

        void createElement(size_t index1,
                           size_t index2,
                           size_t index3,
                           size_t index4,
                           const std::string & materialName);

        void createMoistureBCFixedHc(size_t index1,
                                     size_t index2,
                                     double t_AirTemperature,
                                     double t_ConvectionCoefficient,
                                     double t_Humidity);

        void createMoistureBCVariableHc(size_t index1,
                                        size_t index2,
                                        double t_AirTemperature,
                                        double t_Humidity);

        void createTemperatureBC(size_t index1, size_t index2, double t_Temp1, double t_Temp2);

        void createTemperatureBC(size_t index1, size_t index2, double t_Temp);

        void createBlackBodyRadiationBC(size_t index1,
                                        size_t index2,
                                        double t_Emissivity,
                                        double t_RadiationTemperature);

        void createSimplifiedRadiationBC(size_t index1,
                                         size_t index2,
                                         double t_RadiationCoefficient,
                                         double t_RadiationTemperature);

    private:
        static double normError(const std::vector<double> & vec1, const std::vector<double> & vec2);

        ThermalDomain m_ThermalDomain;
        MoistureDomain m_MoistureDomain;
    };

}   // namespace HygroThermFEM
