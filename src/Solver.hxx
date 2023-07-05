#pragma once

#include <vector>

#include "Element2D.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;
    struct SingleDomain;

    //! \brief Class to hold solution from single timestep.
    struct SingleTimestepSolution
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

    //! Calculates steady state for given data
    //! @param domain Domain for which steady state is being calculated
    //! @param output Domain solution
    std::vector<double> steadyState(SingleDomain & domain);

    //! Calculates steady state solution for multiple domains.
    Solution steadyState(HygroThermFEM::MultiDomain & domain);

    //! \brief Calling timestep calculations
    //! @param previousTimestepStateValues Current state values from previous timestep
    //! @param t_DTime Time different for between timesteps
    //! @param timestepIndex Current timestep index used in variable boundary conditions
    std::pair<std::vector<double>, bool>
      transientTimestep(SingleDomain & domain,
                        const std::vector<double> & previousTimestepStateValues,
                        double t_DTime,
                        size_t timestepIndex);

    //! \brief Calculates next timestep values from current (initial) values
    //! @param domain Domain for which transient solution is being calculated
    //! @param currentStateValues Current values of state variable or initial condition
    //! @param t_DTime Timestep in transient solution
    //! @param timestepIndex Index for current timestep. Used in variable boundary conditions
    //! case. It is defaulted to zero in case of non-variable boundary condition calculations
    //! are requested.
    SingleTimestepSolution transient(SingleDomain & domain,
                                     const std::vector<double> & currentStateValues,
                                     double t_DTime,
                                     size_t timestepIndex = 0);

    //! \brief Calculates next timestep value from current values
    //! \param previousTimestepTemperature vector of nodal temperatures from previous timestep
    //! \param previousTimestepHumidity vector of nodal humidity from previous timestep
    //! \param t_DTime time between two timestep
    //! \param timestepIndex current timestep index used in variable boundary conditions case
    Solution transient(HygroThermFEM::MultiDomain & domain,
                       const std::vector<double> & previousTimestepTemperature,
                       const std::vector<double> & previousTimestepHumidity,
                       double t_DTime,
                       size_t timestepIndex = 0);
}   // namespace HygroThermFEM