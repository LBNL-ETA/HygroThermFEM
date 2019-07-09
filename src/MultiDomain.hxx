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
                 const std::vector<double> & iceContent,
                 const std::vector<NodeFlux> & heatFlux,
                 const std::vector<NodeFlux> & waterFlux,
                 double temperatureError,
                 double humidityError);

        double dTime;
        std::vector<double> temperature;
        std::vector<double> humidity;
        std::vector<double> waterContent;
        std::vector<double> liquidWaterContent;
        std::vector<double> vaporContent;
        std::vector<double> iceContent;
        std::vector<NodeFlux> heatFlux;
        std::vector<NodeFlux> waterFlux;
        double temperatureError;
        double humidityError;
    };

    class MultiDomain
    {
    public:
        MultiDomain(bool performThermal = true, bool performMoisture = true);

        //! \brief Calculates next timestep value from current values
        //! @param temperature vector of nodal temperatures from previous timestep
        //! @param humidity vector of nodal humidity from previous timestep
        //! @param t_DTime time between two timestep
        //! @param timestepIndex current timestep index used in variable boundary conditions case
        Solution transient(std::vector<double> & temperature,
                           std::vector<double> & humidity,
                           double t_DTime,
                           size_t timestepIndex = 0);

        //! Calculates steady state solution for multiple domains.
        Solution steadyState();

        static std::vector<double> property(Variable property);

        //! \brief Creates element with material reference
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        void createElement(size_t index1,
                           size_t index2,
                           size_t index3,
                           size_t index4,
                           const std::string & materialName);

        //! \brief Creates boundary condition with coefficients that are identical during the entire
        //! transient simulation or in steady-state calculations
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBchcCoefficients Structure with values that are used in fixed convection
        //! coefficient boundary condition calculations
        void createMoistureBCFixedHc(size_t index1,
                                     size_t index2,
                                     const FixedBCHCCoefficients & fixedBchcCoefficients);

        //! \brief Creates set of boundary condition coefficients that are used during transient
        //! simulation
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBchcCoefficients Structure with values that are used in fixed convection
        //! coefficient boundary condition calculations 
        void
          createMoistureBCFixedHc(size_t index1,
                                  size_t index2,
                                  const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients);

        void createMoistureBCVariableHc(size_t index1,
                                        size_t index2,
                                        const VariableBCHCCoefficients & varHCCoeff);

        void createTemperatureBC(size_t index1, size_t index2, double t_Temp1, double t_Temp2);

        void createTemperatureBC(size_t index1, size_t index2, double t_Temp);

        void createBlackBodyRadiationBC(size_t index1,
                                        size_t index2,
                                        double t_Emissivity,
                                        double t_RadiationTemperature);

        void createSimplifiedRadiationBC(size_t index1,
                                         size_t index2,
                                         const LinearizedRadiationBCCoefficients & linearRadBC);

    private:
        static double normError(const std::vector<double> & vec1, const std::vector<double> & vec2);

        ThermalDomain m_ThermalDomain;
        MoistureDomain m_MoistureDomain;
        bool m_PerformThermal;
        bool m_PerformMoisture;
    };

}   // namespace HygroThermFEM
