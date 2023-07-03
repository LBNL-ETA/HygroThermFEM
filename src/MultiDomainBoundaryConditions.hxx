#pragma once

#include <vector>

#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;

    //! \brief Creates boundary condition with coefficients that are identical during the entire
    //! transient simulation or in steady-state calculations
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param fixedBchcCoefficients Structure with values that are used in fixed convection
    //! coefficient boundary condition calculations
    void createBC_FixedHc(MultiDomain & domain,
                          size_t index1,
                          size_t index2,
                          const FixedBCHCCoefficients & fixedBchcCoefficients);

    //! \brief Creates set of boundary condition coefficients that are used during transient
    //! simulation
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param fixedBchcCoefficients Structure with values that are used in fixed convection
    //! coefficient boundary condition calculations
    void createBC_FixedHc(MultiDomain & domain,
                          size_t index1,
                          size_t index2,
                          const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients);

    //! \brief Creates boundary condition with coefficients that are identical during the entire
    //! transient simulation or in steady-state calculations
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param varHCCoeff Structure with values that are used in variable convection
    //! coefficient boundary condition calculations
    //! @param Surface tilt at the boundary [degrees]
    void createBC_TARPHc(MultiDomain & domain,
                         size_t index1,
                         size_t index2,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt = 90);

    //! \brief Creates set of boundary condition coefficients that are used during transient
    //! simulation
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param varHCCoeff Structure with values that are used in variable convection
    //! coefficient boundary condition calculations
    //! @param Surface tilt at the boundary [degrees]
    void createBC_TARPHc(MultiDomain & domain,
                         size_t index1,
                         size_t index2,
                         const std::vector<TARPCoefficients> & varHCCoeff,
                         double surfaceTilt = 90);

    //! \brief Creation of ASHRAE inside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    //! @param surfaceHeight Surface height at the boundary [meters]
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    void createBC_ASHRAEInsideHc(MultiDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90);

    //! \brief Creation of ASHRAE inside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    //! @param surfaceHeight Surface height at the boundary [meters]
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    void createBC_ASHRAEInsideHc(MultiDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const std::vector<ASHRAEInsideCoefficients> & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90);

    //! \brief Creation of ASHRAE outside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    void createBC_ASHRAEOutsideHc(MultiDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff);

    //! \brief Creation of ASHRAE outside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    void createBC_ASHRAEOutsideHc(MultiDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<ASHRAEOutsideCoefficients> & coeff);

    //! \brief Creation of Yazdanian-Klems outside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    void createBC_YazdanianKlemsHc(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff);

    //! \brief Creation of Yazdanian-Klems outside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    void createBC_YazdanianKlemsHc(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<YazdanianKlemsCoefficients> & coeff);

    //! \brief Creation of Kimura outside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    void createBC_KimuraHc(MultiDomain & domain,
                           size_t index1,
                           size_t index2,
                           const std::vector<KimuraCoefficients> & coeff);

    //! \brief Creation of Kimura outside convection boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    void createBC_KimuraHc(MultiDomain & domain,
                           size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff);

    //! \brief Sets fixed temperature boundary conditions
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param temp1 Temperature value that will be set at the node 1
    //! @param temp2 Temperature value that will be set at the node 2
    void createBC_FixedTemperature(
      MultiDomain & domain, size_t index1, size_t index2, double temp1, double temp2);

    //! \brief Sets fixed temperature boundary conditions
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param temp Set of node temperatures for every given timstep (each node can have
    //! different temperature).
    void createBC_FixedTemperature(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<ConstantBCTemperatures> & temp);

    //! \brief Sets fixed temperature boundary conditions
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param temp Temperature value that will be set at the nodes
    void createBC_FixedTemperature(MultiDomain & domain, size_t index1, size_t index2, double temp);

    //! \brief Sets fixed temperature boundary conditions at every timestep
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param temp Temperature values at every timestep that will be set at the nodes
    void createBC_FixedTemperature(MultiDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   std::vector<double> temp);

    //! \brief Sets constant temperature and humidity at the boundary
    //!
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param values Temperature and humidity values that will be set at the nodes
    void createBC_FixedTemperatureAndHumidity(MultiDomain & domain,
                                              size_t index1,
                                              size_t index2,
                                              const TemperatureAndHumidity & values);

    //! \brief Sets constant temperature and humidity at the boundary
    //!
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param values Temperature and humidity values that will be set at the nodes
    void createBC_FixedTemperatureAndHumidity(MultiDomain & domain,
                                              size_t index1,
                                              size_t index2,
                                              const std::vector<TemperatureAndHumidity> & values);

    //! \brief Creation of flux boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Flux Constant flux in both nodes
    void createBC_FixedHeatFlux(MultiDomain & domain, size_t index1, size_t index2, double t_Flux);

    //! \brief Creation of flux boundary condition
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Flux Constant flux in both nodes
    void createBC_FixedHeatFlux(MultiDomain & domain,
                                size_t index1,
                                size_t index2,
                                const std::vector<double> & t_Flux);

    //! \brief Sets radiation boundary condition that is fixed during entire transient
    //! simulation or set for steady-state case
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param t_Emissivity Surface emissivity
    //! @param t_RadiationTemperature Environment radiation temperature
    void createBC_BlackBodyRadiation(MultiDomain & domain,
                                     size_t index1,
                                     size_t index2,
                                     double t_Emissivity,
                                     double t_RadiationTemperature);

    //! \brief Sets radiation boundary condition that is fixed during entire transient
    //! simulation or set for steady-state case
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param radCoeffs Radiation coefficients for every timestep
    void
      createBC_BlackBodyRadiation(MultiDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs);

    //! \brief Sets radiation boundary condition that is fixed during entire transient
    //! simulation or set for steady-state case
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param radCoeffs Radiation coefficients for entire transient simulation or steady-state
    void createBC_LinearizedRadiation(MultiDomain & domain,
                                      size_t index1,
                                      size_t index2,
                                      const LinearizedRadiationBCCoefficients & linearRadBC);

    //! \brief Sets radiation boundary condition that is fixed during entire transient
    //! simulation or set for steady-state case
    //! @param domain Domain at which boundary conditions will be applied
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param radCoeffs Radiation coefficients for every timestep
    void createBC_LinearizedRadiation(
      MultiDomain & domain,
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC);
}   // namespace HygroThermFEM