#pragma once

#include <utility>
#include <vector>

#include "Element2D.hxx"

namespace HygroThermFEM
{
    //! \brief Class to hold solution from single timestep.
    struct SingleTimestepSolution
    {
        SingleTimestepSolution(const std::vector<double> & solution, double dtime) :
            dTime(dtime),
            solution(solution)
        {}
            std::vector<double> solution;   //!< Solution from single transient step
            double dTime;   //!< Timestep for which solution has been performed. Engine can adopt
                            //!< new timestep for which solution will converge.
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
    }   // namespace HygroThermFEM