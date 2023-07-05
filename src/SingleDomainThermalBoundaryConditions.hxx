#pragma once

#include <vector>
#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    struct SingleDomain;
}

namespace HygroThermFEM::Thermal
{
    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param fixedBCHCCoefficients structure to hold fixed convection coefficient boundary
    //! conditions
    //! @param t_CalculateMoisture Flag on whether or not to calculate moisture
    void createBC_FixedHc(HygroThermFEM::SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const FixedBCHCCoefficients & fixedBCHCCoefficients,
                          bool t_CalculateMoisture = true);

    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param fixedBCHCCoefficients structure to hold fixed convection coefficients boundary
    //! conditions for transient simulations
    //! @param calculateMoisture Flag on whether or not to calculate moisture
    void createBC_FixedHc(HygroThermFEM::SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const std::vector<FixedBCHCCoefficients> & fixedBCHCCoefficients,
                          bool calculateMoisture = true);

    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param varHCCoeff structure to hold variable convection coefficient boundary conditions
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    //! @param simulateVaporFluxEnergy Flag on whether or not to include energy from moisture
    //! flux
    void createBC_TARPHc(HygroThermFEM::SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt = 90,
                         bool simulateVaporFluxEnergy = true);

    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param varHCCoeff structure to hold variable convection coefficients boundary
    //! conditions for transient simulations
    //! @param surfaceTilt Surface tilt at the current boundary
    //! @param simulateVaporFluxEnergy Flag on whether or not to calculate moisture
    void createBC_TARPHc(HygroThermFEM::SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const std::vector<TARPCoefficients> & varHCCoeff,
                         double surfaceTilt = 90,
                         bool simulateVaporFluxEnergy = true);

    //! \brief Creation of ASHRAE inside convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    //! @param surfaceHeight Surface height at the boundary [meters]
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    //! @param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_ASHRAEInsideHc(HygroThermFEM::SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90,
                                 bool simulateVaporFluxEnergy = true);

    //! \brief Creation of ASHRAE inside convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    //! @param surfaceHeight Surface height at the boundary [meters]
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    //! @param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_ASHRAEInsideHc(HygroThermFEM::SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const std::vector<ASHRAEInsideCoefficients> & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90,
                                 bool simulateVaporFluxEnergy = true);

    //! \brief Creation of ASHRAE outside convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable through timesteps
    //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_ASHRAEOutsideHc(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff,
                                  bool simulateVaporFluxEnergy = true);

    //! \brief Creation of ASHRAE outside convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable through timesteps
    //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_ASHRAEOutsideHc(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<ASHRAEOutsideCoefficients> & coeff,
                                  bool simulateVaporFluxEnergy = true);

    //! \brief Creation of Yazdanian-Klems convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable through timesteps
    //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_YazdanianKlemsHc(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff,
                                   bool simulateVaporFluxEnergy = true);

    //! \brief Creation of Yazdanian-Klems convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable through timesteps
    //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_YazdanianKlemsHc(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<YazdanianKlemsCoefficients> & coeff,
                                   bool simulateVaporFluxEnergy = true);

    //! \brief Creation of Kimura convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable through timesteps
    //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_KimuraHc(HygroThermFEM::SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff,
                           bool simulateVaporFluxEnergy = true);

    //! \brief Creation of Kimura convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable through timesteps
    //! \param simulateVaporFluxEnergy Flag to indicate whether or not energy from vapor flux
    //! will be calculated.
    void createBC_KimuraHc(HygroThermFEM::SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const std::vector<KimuraCoefficients> & coeff,
                           bool simulateVaporFluxEnergy = true);

    //! \brief Creation of temperature boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Temp1 Constant temperature at node 1
    //! @param t_Temp2 Constant temperature at node 2
    void createBC_FixedTemperature(
      HygroThermFEM::SingleDomain & domain, size_t index1, size_t index2, double t_Temp1, double t_Temp2);


    //! \brief Sets fixed temperature boundary conditions
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param temp Set of node temperatures for every given timstep (each node can have
    //! different temperature).
    void createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<ConstantBCTemperatures> & temp);

    //! \brief Creation of temperature boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Temp Constant temperature in both nodes
    void
      createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain, size_t index1, size_t index2, double t_Temp);

    //! \brief Creation of temperature boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Temp Constant temperature at each timestep in both nodes
    void createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   std::vector<double> t_Temp);

    //! \brief Creation of flux boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Flux Constant flux in both nodes
    void createBC_FixedFlux(HygroThermFEM::SingleDomain & domain, size_t index1, size_t index2, double t_Flux);

    //! \brief Creation of flux boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Flux Constant flux in both nodes at every timestep
    void createBC_FixedFlux(HygroThermFEM::SingleDomain & domain,
                            size_t index1,
                            size_t index2,
                            const std::vector<double> & t_Flux);

    //! \brief Creation of black body radiation
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Emissivity Emissivity at boundary condition
    //! @param t_RadiationTemperature Radiation temperature
    void createBC_BlackBodyRadiation(HygroThermFEM::SingleDomain & domain,
                                     size_t index1,
                                     size_t index2,
                                     double t_Emissivity,
                                     double t_RadiationTemperature);

    //! \brief Sets radiation boundary condition that is fixed during entire transient
    //! simulation or set for steady-state case
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param radCoeffs Radiation coefficients for every timestep
    void
      createBC_BlackBodyRadiation(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs);

    //! \brief Creation of linearized radiation boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param linearRadBC Structure that holds coefficients for linearized radiation
    void createBC_LinearizedRadiation(HygroThermFEM::SingleDomain & domain,
                                      size_t index1,
                                      size_t index2,
                                      const LinearizedRadiationBCCoefficients & linearRadBC);

    //! \brief Creation of linearized radiation boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param linearRadBC Linarized radiation boundary condition at every timestep
    void createBC_LinearizedRadiation(
      HygroThermFEM::SingleDomain & domain,
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC);
}   // namespace HygroThermFEM