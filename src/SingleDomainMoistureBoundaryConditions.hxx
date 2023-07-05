#pragma once

#include <vector>
#include "BoundaryConditionCoefficients.hxx"

namespace HygroThermFEM
{
    struct SingleDomain;
}

namespace HygroThermFEM::Moisture
{
    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param varHCCoeff structure to hold variable convection coefficient boundary conditions
    //! @param surfaceTilt Surface tilt at the boundary
    void createBC_TARPHc(HygroThermFEM::SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt = 90);

    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param varCoeff structure to hold fixed convection coefficient boundary
    //! conditions for every timestep
    //! @param surfaceTilt Surface tilt at the boundary
    void createBC_TARPHc(HygroThermFEM::SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const std::vector<TARPCoefficients> & varCoeff,
                         double surfaceTilt = 90);

    //! \brief Creation of ASHRAE inside convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    //! @param surfaceHeight Surface height at the boundary [meters]
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    void createBC_ASHRAEInsideHc(HygroThermFEM::SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90);

    //! \brief Creation of ASHRAE inside convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param coeff Structure to hold variable convection coefficients that are variable
    //! through timesteps
    //! @param surfaceHeight Surface height at the boundary [meters]
    //! @param surfaceTilt Surface tilt at the boundary [degrees]
    void createBC_ASHRAEInsideHc(HygroThermFEM::SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const std::vector<ASHRAEInsideCoefficients> & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt = 90);

    //! \brief Creation of ASHRAE outside convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable at different timesteps
    void createBC_ASHRAEOutsideHc(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff);

    //! \brief Creation of ASHRAE outside convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable at different timesteps
    void createBC_ASHRAEOutsideHc(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<ASHRAEOutsideCoefficients> & coeff);

    //! \brief Creation of Yazdanian-Klems convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable at different timesteps
    void createBC_YazdanianKlemsHc(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff);

    //! \brief Creation of Yazdanian-Klems convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable at different timesteps
    void createBC_YazdanianKlemsHc(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<YazdanianKlemsCoefficients> & coeff);

    //! \brief Creation of Kimura convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable at different timesteps
    void createBC_KimuraHc(HygroThermFEM::SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff);


    //! \brief Creation of Kimura convection boundary condition
    //!
    //! @param domain Domain to add boundary condition to
    //! \param index1 Node 1 index
    //! \param index2 Node 2 index
    //! \param coeff Structure to hold coefficients that are variable at different timesteps
    void createBC_KimuraHc(HygroThermFEM::SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const std::vector<KimuraCoefficients> & coeff);

    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param fixedBchcCoefficients structure to hold fixed convection coefficient boundary
    //! conditions
    void createBC_FixedHc(HygroThermFEM::SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const FixedBCHCCoefficients & fixedBchcCoefficients);

    //! \brief Creation of convection boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param fixedBchcCoefficients structure to hold fixed convection coefficient boundary
    //! conditions for every timestep
    void createBC_FixedHc(HygroThermFEM::SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients);

    //! \brief Creation constant temperature and humidity boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param values structure to hold fixed temperature and humidity
    void createBC_FixedHumidity(HygroThermFEM::SingleDomain & domain,
                                size_t index1,
                                size_t index2,
                                const TemperatureAndHumidity & values);

    //! \brief Creation constant temperature and humidity boundary condition
    //! @param domain Domain to add boundary condition to
    //! @param index1 Node 1 index
    //! @param index2 Node 2 index
    //! @param values structure to hold fixed temperature and humidity for every timestep
    void createBC_FixedHumidity(HygroThermFEM::SingleDomain & domain,
                                size_t index1,
                                size_t index2,
                                const std::vector<TemperatureAndHumidity> & values);
}   // namespace HygroThermFEM::Moisture