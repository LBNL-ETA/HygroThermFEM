#pragma once

#include <vector>

#include "Element2D.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;

    //! \brief Class to hold solution from single timestep.
    struct SingleSolution
    {
        std::vector<double> solution;   //!< Solution from single transient step
        double dTime;   //!< Timestep for which solution has been performed. Engine can adopt new
                                        //!< timestep for which solution will converge.
    };

    //! \brief Keeps solution from current timestep for every node in the domain.
    struct Solution
    {
        Solution(double dtime,
                 std::vector<double> temperature,
                 std::vector<double> humidity,
                 std::vector<double> waterContent,
                 std::vector<double> liquidWaterContent,
                 std::vector<double> vaporContent,
                 std::vector<double> iceContent,
                 std::vector<NodeFlux> heatFlux,
                 std::vector<NodeFlux> waterFlux,
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

    //! \brief Calculates next timestep value from current values
    //! \param previousTimestepTemperature vector of nodal temperatures from previous timestep
    //! \param previousTimestepHumidity vector of nodal humidity from previous timestep
    //! \param t_DTime time between two timestep
    //! \param timestepIndex current timestep index used in variable boundary conditions case
    Solution transient(HygroThermFEM::MultiDomain & domain, const std::vector<double> & previousTimestepTemperature,
                       const std::vector<double> & previousTimestepHumidity,
                       double t_DTime,
                       size_t timestepIndex = 0);

    //! Calculates steady state solution for multiple domains.
    Solution steadyState(HygroThermFEM::MultiDomain & domain);
}   // namespace HygroThermFEM