#pragma once

#include <memory>
#include <vector>

#include "SquareMatrix.hxx"
#include "IBCLine2D.hxx"

namespace HygroThermFEM
{
    //! \brief Container for all boundary conditions.
    //!
    //! Beside simple container that will hold all kinds of boundary conditions, responsibility of
    //! this class is to create assembled vector and matrix for all boundary conditions.
    //! Boundary conditions can also be inserted for transient case where each timestep will have
    //! specific boundary condition.
    class BoundaryConditions2D
    {
    public:
        //! Simple constructor
        explicit BoundaryConditions2D();

        //! Making copy of boundary conditions container is not allowed. Since all changes
        //! and updates should be kept at one place.
        BoundaryConditions2D(const BoundaryConditions2D & other) = delete;

        //! \brief Assembled H matrix from all boundary conditions.
        //! @param timestepIndex Timestep index for which matrix being calculated. It is defaulted
        //! to zero for single boundary condition input
        SquareMatrix HMatrix(size_t timestepIndex = 0) const;

        //! \brief Assembled R vector from all boundary conditions.
        //! @param timestepIndex Timestep index for which matrix being calculated. It is defaulted
        //! to zero for single boundary condition input
        std::vector<double> RVector(size_t timestepIndex = 0) const;

        //! Returns linearity of the problem. It is important for domain to know what set of
        //! equations should be applied.
        bool isLinear() const;

        //! Assigns new steady-state or single transient boundary condition to the pool
        //! @param bc single boundary condition used for steady-state or transient
        void assignBC(std::unique_ptr<IBCLinear2D> && bc);

        //! Assigns new transient boundary condition for next timestep. Timestep index is not
        //! traced. When simulation starts, it is simple requirement that number of boundary
        //! condition must match number of timesteps.
        //! @param bc multiple timestep boundary conditions (different values for each timestep)
        void assignTimestepBCs(std::vector<std::unique_ptr<IBCLinear2D>> && bc);

    protected:
        // Boundary Conditions used either for steady-state or transient cases.
        std::vector<std::unique_ptr<IBCLinear2D>> m_BCs;

        // Boundary conditions that are used only in transient cases where inside vector means
        // different boundary condition for every timestep.
        std::vector<std::vector<std::unique_ptr<IBCLinear2D>>> m_TransientBCs;
        bool m_Linear;
    };

}   // namespace HygroThermFEM
