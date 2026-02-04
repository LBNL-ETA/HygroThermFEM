#pragma once

#include "Domain.hxx"

namespace HygroThermFEM
{
    class MultiDomain;

    //! \brief Domain class for solving temperature solution.
    class ThermalDomain : public IDomain
    {
        friend class MultiDomain;

    private:
        //! Simple constructor - only accessible via MultiDomain
        explicit ThermalDomain(Nodes & nodePool,
                               Materials & materialPool,
                               bool automaticUpdatePreviousTimestep = true);

        //! \brief Creates and adds element into domain - only accessible via MultiDomain
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        void createElement(size_t index1,
                           size_t index2,
                           size_t index3,
                           size_t index4,
                           const std::string & materialName) override;

    public:

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBCHCCoefficients structure to hold fixed convection coefficient boundary
        //! conditions
        //! @param t_CalculateMoisture Flag on whether or not to calculate moisture
        void createBC_FixedHc(size_t index1,
                              size_t index2,
                              const FixedBCHCCoefficients & fixedBCHCCoefficients,
                              bool t_CalculateMoisture = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBCHCCoefficients structure to hold fixed convection coefficients boundary
        //! conditions for transient simulations
        //! @param calculateMoisture Flag on whether or not to calculate moisture
        void createBC_FixedHc(size_t index1,
                              size_t index2,
                              const std::vector<FixedBCHCCoefficients> & fixedBCHCCoefficients,
                              bool calculateMoisture = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff structure to hold variable convection coefficient boundary conditions
        //! @param surfaceTilt Surface tilt at the boundary [degrees]
        //! @param simulateVaporFluxEnergy Flag on whether or not to include energy from moisture
        //! flux
        void createBC_TARPHc(size_t index1,
                             size_t index2,
                             const TARPCoefficients & varHCCoeff,
                             double surfaceTilt = 90,
                             bool simulateVaporFluxEnergy = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff structure to hold variable convection coefficients boundary
        //! conditions for transient simulations
        //! @param surfaceTilt Surface tilt at the current boundary
        //! @param simulateVaporFluxEnergy Flag on whether or not to calculate moisture
        void createBC_TARPHc(size_t index1,
                             size_t index2,
                             const std::vector<TARPCoefficients> & varHCCoeff,
                             double surfaceTilt = 90,
                             bool simulateVaporFluxEnergy = true);

        //! \brief Creation of ASHRAE inside convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param coeff Structure to hold variable convection coefficients that are variable
        //! through timesteps
        //! @param surfaceHeight Surface height at the boundary [meters]
        //! @param surfaceTilt Surface tilt at the boundary [degrees]
        //! @param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_ASHRAEInsideHc(size_t index1,
                                     size_t index2,
                                     const ASHRAEInsideCoefficients & coeff,
                                     double surfaceHeight,
                                     double surfaceTilt = 90,
                                     bool simulateVaporFluxEnergy = true);

        //! \brief Creation of ASHRAE inside convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param coeff Structure to hold variable convection coefficients that are variable
        //! through timesteps
        //! @param surfaceHeight Surface height at the boundary [meters]
        //! @param surfaceTilt Surface tilt at the boundary [degrees]
        //! @param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_ASHRAEInsideHc(size_t index1,
                                     size_t index2,
                                     const std::vector<ASHRAEInsideCoefficients> & coeff,
                                     double surfaceHeight,
                                     double surfaceTilt = 90,
                                     bool simulateVaporFluxEnergy = true);

        //! \brief Creation of ASHRAE outside convection boundary condition
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Structure to hold coefficients that are variable through timesteps
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_ASHRAEOutsideHc(size_t index1,
                                      size_t index2,
                                      const ASHRAEOutsideCoefficients & coeff,
                                      bool simulateVaporFluxEnergy = true);

        //! \brief Creation of ASHRAE outside convection boundary condition
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Structure to hold coefficients that are variable through timesteps
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_ASHRAEOutsideHc(size_t index1,
                                      size_t index2,
                                      const std::vector<ASHRAEOutsideCoefficients> & coeff,
                                      bool simulateVaporFluxEnergy = true);

        //! \brief Creation of Yazdanian-Klems convection boundary condition
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Structure to hold coefficients that are variable through timesteps
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_YazdanianKlemsHc(size_t index1,
                                      size_t index2,
                                      const YazdanianKlemsCoefficients & coeff,
                                      bool simulateVaporFluxEnergy = true);

        //! \brief Creation of Yazdanian-Klems convection boundary condition
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Structure to hold coefficients that are variable through timesteps
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_YazdanianKlemsHc(size_t index1,
                                      size_t index2,
                                      const std::vector<YazdanianKlemsCoefficients> & coeff,
                                      bool simulateVaporFluxEnergy = true);

        //! \brief Creation of Kimura convection boundary condition
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Structure to hold coefficients that are variable through timesteps
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_KimuraHc(size_t index1,
                                       size_t index2,
                                       const KimuraCoefficients & coeff,
                                       bool simulateVaporFluxEnergy = true);

        //! \brief Creation of Kimura convection boundary condition
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param coeff Structure to hold coefficients that are variable through timesteps
        //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
        //! will be calculated.
        void createBC_KimuraHc(size_t index1,
                               size_t index2,
                               const std::vector<KimuraCoefficients> & coeff,
                               bool simulateVaporFluxEnergy = true);

        //! \brief Creation of temperature boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Temp1 Constant temperature at node 1
        //! @param t_Temp2 Constant temperature at node 2
        void
          createBC_FixedTemperature(size_t index1, size_t index2, double t_Temp1, double t_Temp2);


        //! \brief Sets fixed temperature boundary conditions
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param temp Set of node temperatures for every given timstep (each node can have
        //! different temperature).
        void createBC_FixedTemperature(size_t index1,
                                       size_t index2,
                                       const std::vector<ConstantBCTemperatures> & temp);

        //! \brief Creation of temperature boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Temp Constant temperature in both nodes
        void createBC_FixedTemperature(size_t index1, size_t index2, double t_Temp);

        //! \brief Creation of temperature boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Temp Constant temperature at each timestep in both nodes
        void createBC_FixedTemperature(size_t index1, size_t index2, std::vector<double> t_Temp);

        //! \brief Creation of flux boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Flux Constant flux in both nodes
        void createBC_FixedFlux(size_t index1, size_t index2, double t_Flux);

        //! \brief Creation of flux boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Flux Constant flux in both nodes at every timestep
        void createBC_FixedFlux(size_t index1, size_t index2, std::vector<double> t_Flux);

        //! \brief Creation of black body radiation
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param t_Emissivity Emissivity at boundary condition
        //! @param t_RadiationTemperature Radiation temperature
        void createBC_BlackBodyRadiation(size_t index1,
                                         size_t index2,
                                         double t_Emissivity,
                                         double t_RadiationTemperature);

        //! \brief Sets radiation boundary condition that is fixed during entire transient
        //! simulation or set for steady-state case
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param radCoeffs Radiation coefficients for every timestep
        void createBC_BlackBodyRadiation(
          size_t index1,
          size_t index2,
          const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs);

        //! \brief Creation of linearized radiation boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param linearRadBC Structure that holds coefficients for linearized radiation
        void createBC_LinearizedRadiation(size_t index1,
                                          size_t index2,
                                          const LinearizedRadiationBCCoefficients & linearRadBC);

        //! \brief Creation of linearized radiation boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param linearRadBC Linarized radiation boundary condition at every timestep
        void createBC_LinearizedRadiation(
          size_t index1,
          size_t index2,
          const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC);

    protected:
        void postProcess(std::vector<double> & solution) override;
    };
}   // namespace HygroThermFEM
