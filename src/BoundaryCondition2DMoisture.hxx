#pragma once
#include "BoundaryCondition2D.hxx"

namespace HygroThermFEM
{
    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient with TARP
    //! algorithm.
    class MoistureBCTARPHc : public IMoistureBC
    {
    public:
        //! \brief Construction of boundary condition with TARP algorithm for convective heat
        //! transfer calculations
        //!
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param materialName Material name for which boundary is associated with
        //! @param varHCCoeff Coefficients that are necessary for TARP heat transfer coefficient
        //! calculations. Structure contain only coefficients that are variable through every
        //! timestep.
        //! @param surafceTilt Surface tilt at the boundary. [degrees]
        MoistureBCTARPHc(size_t index1,
                         size_t index2,
                         const std::string & materialName,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt = 90);
    };

    /////////////////////////////////////////////////////
    /// MoistureBCASHRAEInside
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient based on ASHRAE
    //! inside algorithm during one timestep or steady-state.
    class MoistureBCASHRAEInside : public IMoistureBC
    {
    public:
        //! \brief Moisture boundary condition that calculates convective coefficient based on
        //! ASHRAE inside algorithm
        //!
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param materialName Material that is adjacent to the boundary
        //! @param coeffs Coefficients that are necessary for ASHRAE inside convection heat transfer
        //! coefficient calculation Structure only contain coefficients that are variable through
        //! every timestep.
        MoistureBCASHRAEInside(size_t index1,
                               size_t index2,
                               const std::string & materialName,
                               const ASHRAEInsideCoefficients & coeffs,
                               double surfaceHeight,
                               double surfaceTilt = 90);
    };

    /////////////////////////////////////////////////////
    /// MoistureBCASHRAEOutside
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient based on ASHRAE
    //! outside algorithm during one timestep or steady-state.
    class MoistureBCASHRAEOutside : public IMoistureBC
    {
    public:
        //! \brief Moisture boundary condition that calculates convective coefficient based on
        //! ASHRAE outside algorithm.
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param materialName Material that is adjacent to the boundary
        //! \param coeffs Coefficients that are necessary for ASHRAE outside convection heat
        //! transfer coefficients calculation. Structure only contain coefficients that are variable
        //! through timestep.
        MoistureBCASHRAEOutside(size_t index1,
                                size_t index2,
                                const std::string & materialName,
                                const ASHRAEOutsideCoefficients & coeffs);
    };

    /////////////////////////////////////////////////////
    /// YazdanianKlemsBC
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient based on ASHRAE
    //! outside algorithm during one timestep or steady-state.
    class MoistureYazdanianKlemsBC : public IMoistureBC
    {
    public:
        //! \brief Moisture boundary condition that calculates convective coefficient based on
        //! Yazdanian-Klems outside algorithm.
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param materialName Material that is adjacent to the boundary
        //! \param coeffs Coefficients that are necessary for Yazdanian-Klems outside convection
        //! heat transfer coefficients calculation. Structure only contain coefficients that are
        //! variable through timestep.
        MoistureYazdanianKlemsBC(size_t index1,
                                 size_t index2,
                                 const std::string & materialName,
                                 const YazdanianKlemsCoefficients & coeffs);
    };

    //! \brief Moisture boundary condition that calculates convective coefficient based on Kimura
    //! algorithm during one timestep or steady-state.
    class MoistureKimuraBC : public IMoistureBC
    {
    public:
        //! \brief Moisture boundary condition that calculates convective coefficient based on
        //! Kimura model.
        //!
        //! \param index1 Node 1 index
        //! \param index2 Node 2 index
        //! \param materialName Material that is adjacent to the boundary
        //! \param coeffs Coefficients that are necessary for Kmura outside convection heat transfer
        //! coefficient calculation. Structure only contain coefficients that are variable through
        //! timestep.
        MoistureKimuraBC(size_t index1,
                         size_t index2,
                         const std::string & materialName,
                         const KimuraCoefficients & coeffs);
    };

    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition with fixed convective coefficient
    class MoistureBCFixedHc : public IMoistureBC
    {
    public:
        MoistureBCFixedHc(size_t index1,
                          size_t index2,
                          const std::string & materialName,
                          const FixedBCHCCoefficients & fixedBchcCoefficients);
    };

    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////

    class MoistureBCFixedHumidity : public MoistureBCFixedHc
    {
    public:
        MoistureBCFixedHumidity(size_t index1,
                                size_t index2,
                                const std::string & materialName,
                                TemperatureAndHumidity values) :
            MoistureBCFixedHc(index1,
                              index2,
                              materialName,
                              {values.Temperature, hugeFilmCoefficient, values.Humidity})
        {}

    private:
        inline static const double hugeFilmCoefficient{1e18};
    };

}   // namespace HygroThermFEM
