#pragma once

#include "Domain.hxx"

namespace HygroThermFEM
{
    //! \brief Domain class for solving humidity distribution.
    class MoistureDomain : public IDomain
    {
    public:
        //! Simple constructor
        MoistureDomain(bool automaticUpdatePreviousTimestep = true);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varHCCoeff structure to hold variable convection coefficient boundary conditions
        //! @param surfaceTilt Surface tilt at the boundary
        void createBC_TARPHc(size_t index1,
                             size_t index2,
                             const TARPCoefficients & varHCCoeff,
                             double surfaceTilt = 90);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param varCoeff structure to hold fixed convection coefficient boundary
        //! conditions for every timestep
        //! @param surfaceTilt Surface tilt at the boundary
        void createBC_TARPHc(size_t index1,
                             size_t index2,
                             const std::vector<TARPCoefficients> & varCoeff,
                             double surfaceTilt = 90);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBchcCoefficients structure to hold fixed convection coefficient boundary
        //! conditions
        void createBC_FixedHc(size_t index1,
                              size_t index2,
                              const FixedBCHCCoefficients & fixedBchcCoefficients);

        //! \brief Creation of convection boundary condition
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param fixedBchcCoefficients structure to hold fixed convection coefficient boundary
        //! conditions for every timestep
        void createBC_FixedHc(size_t index1,
                              size_t index2,
                              const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients);

        //! \brief Creates and adds element into domain.
        //! @param index1 Node 1 index
        //! @param index2 Node 2 index
        //! @param index3 Node 3 index
        //! @param index4 Node 4 index
        //! @param materialName Material name assigned to the element
        virtual void createElement(size_t index1,
                                   size_t index2,
                                   size_t index3,
                                   size_t index4,
                                   const std::string & materialName) override;

    protected:
        void postProcess(std::vector<double> & solution) override;
    };
}
