#pragma once

#include "Domains.hxx"
#include "MaterialMissingProperties.hxx"

namespace Timesteps
{
    class TimestepObserver;
}

namespace HygroThermFEM
{
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

    class MultiDomain
    {
    public:
        MultiDomain(bool performThermal = true, bool performMoisture = true);

        //! \brief Calculates next timestep value from current values
        //! \param temperature vector of nodal temperatures from previous timestep
        //! \param humidity vector of nodal humidity from previous timestep
        //! \param t_DTime time between two timestep
        //! \param timestepIndex current timestep index used in variable boundary conditions case
        Solution transient(std::vector<double> & temperature,
                           std::vector<double> & humidity,
                           double t_DTime,
                           size_t timestepIndex = 0);

        //! Calculates steady state solution for multiple domains.
        Solution steadyState();

        //! \brief Checks validity of materials for any simulation
        //!
        //! \param simulationType Enumerator for simulation type (SteadyState or Transient)
        //! \return Information about materials with all missing properties.
        MaterialsErrorCheckVector checkForMaterialsValidity(SimulationType simulationType) const;

        //! \brief Sets if moisture simulation will be performed
        //!
        //! @param val if set to true, moisture simulation will be performed
        void performMoistureSimulation(bool val);

        //! \brief Checks if moisture simulation is ON
        bool isMoistureSimulationON() const;

        //! \brief Sets if thermal simulation will be performed
        //!
        //! @param val if set to true, thermal simulation will be performed
        void performThermalSimulation(bool val);

        //! \brief Checks if thermal simulation is ON
        bool isThermalSimulationON() const;

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

        //! \brief Creates boundary condition with coefficients that are identical during the entire
        //! transient simulation or in steady-state calculations
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff Structure with values that are used in variable convection
        //! coefficient boundary condition calculations
        void createMoistureBCVariableHc(size_t index1,
                                        size_t index2,
                                        const VariableBCHCCoefficients & varHCCoeff);

        //! \brief Creates set of boundary condition coefficients that are used during transient
        //! simulation
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff Structure with values that are used in variable convection
        //! coefficient boundary condition calculations
        void createMoistureBCVariableHc(size_t index1,
                                        size_t index2,
                                        const std::vector<VariableBCHCCoefficients> & varHCCoeff);

        //! \brief Sets fixed temperature boundary conditions
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param temp1 Temperature value that will be set at the node 1
        //! @param temp2 Temperature value that will be set at the node 2
        void createTemperatureBC(size_t index1, size_t index2, double temp1, double temp2);

        //! \brief Sets fixed temperature boundary conditions
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param temp Set of node temperatures for every given timstep (each node can have
        //! different temperature).
        void createTemperatureBC(size_t index1,
                                 size_t index2,
                                 const std::vector<ConstantBCTemperatures> & temp);

        //! \brief Sets fixed temperature boundary conditions
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param temp Temperature value that will be set at the nodes
        void createTemperatureBC(size_t index1, size_t index2, double temp);

        //! \brief Sets fixed temperature boundary conditions at every timestep
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param temp Temperature values at every timestep that will be set at the nodes
        void createTemperatureBC(size_t index1, size_t index2, std::vector<double> temp);

        //! \brief Sets radiation boundary condition that is fixed during entire transient
        //! simulation or set for steady-state case
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Emissivity Surface emissivity
        //! @param t_RadiationTemperature Environment radiation temperature
        void createBlackBodyRadiationBC(size_t index1,
                                        size_t index2,
                                        double t_Emissivity,
                                        double t_RadiationTemperature);

        //! \brief Sets radiation boundary condition that is fixed during entire transient
        //! simulation or set for steady-state case
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param radCoeffs Radiation coefficients for every timestep
        void createBlackBodyRadiationBC(
          size_t index1,
          size_t index2,
          const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs);

        //! \brief Sets radiation boundary condition that is fixed during entire transient
        //! simulation or set for steady-state case
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param radCoeffs Radiation coefficients for entire transient simulation or steady-state
        void createLinearizedRadiationBC(size_t index1,
                                         size_t index2,
                                         const LinearizedRadiationBCCoefficients & linearRadBC);

        //! \brief Sets radiation boundary condition that is fixed during entire transient
        //! simulation or set for steady-state case
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param radCoeffs Radiation coefficients for every timestep
        void createLinearizedRadiationBC(
          size_t index1,
          size_t index2,
          const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC);

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! @param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

        //! \brief Assign observer to thermal part of the engine.
        //!
        //! \param observer
        void subscribeThermal(Timesteps::TimestepObserver * observer);

        //! \brief Unsubscribe from thermal notifications
        //!
        //! \param observer
        void unsubscribeThermal(Timesteps::TimestepObserver * observer);

        //! \brief Assign observer to moisture part of the engine.
        //!
        //! \param observer
        void subscribeMoisture(Timesteps::TimestepObserver * observer);

        //! \brief Unsubscribe from moisture notifications
        //!
        //! \param observer
        void unsubscribeMoisture(Timesteps::TimestepObserver * observer);


    private:
        //! \brief Checks validity of materials for transient simulation
        MaterialsErrorCheckVector checkMaterialsForTransientSimulation() const;

        //! \brief Checks validity of materials for transient simulation
        MaterialsErrorCheckVector checkMaterialsForSteadyStateSimulation() const;

        static double normError(const std::vector<double> & vec1, const std::vector<double> & vec2);

        ThermalDomain m_ThermalDomain;
        MoistureDomain m_MoistureDomain;
        bool m_PerformThermal;
        bool m_PerformMoisture;
    };

}   // namespace HygroThermFEM
