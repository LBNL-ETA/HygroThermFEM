#pragma once
#include "BoundaryCondition2D.hxx"

namespace HygroThermFEM
{
    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////

    //! \brief Moisture boundary condition that calculates convective coefficient.
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
}
