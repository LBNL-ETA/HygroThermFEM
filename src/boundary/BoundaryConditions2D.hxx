#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

#include "SquareMatrix.hxx"
#include "IBCLine2D.hxx"

namespace HygroThermFEM
{
    //! \brief Heat rate through one boundary-condition element (a two-node boundary segment).
    //! Positive when heat leaves the domain through the segment (outward-normal convention).
    struct BoundaryHeatRate
    {
        size_t node1{0u};
        size_t node2{0u};
        double heatRate{0.0};
    };

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

        //! Clearing boundary conditions
        void clear();

        //! Making copy of boundary conditions container is not allowed. Since all changes
        //! and updates should be kept at one place.
        BoundaryConditions2D(const BoundaryConditions2D & other) = delete;

        //! \brief Assembled H matrix from all boundary conditions.
        //! @param maxNodeIndex Maximum node index for matrix sizing.
        //! @param timestepIndex Timestep index for which matrix being calculated. It is defaulted
        //! to zero for single boundary condition input
        SquareMatrix HMatrix(size_t maxNodeIndex, size_t timestepIndex = 0) const;

        //! \brief Assembled R vector from all boundary conditions.
        //! @param maxNodeIndex Maximum node index for vector sizing.
        //! @param timestepIndex Timestep index for which matrix being calculated. It is defaulted
        //! to zero for single boundary condition input
        std::vector<double> RVector(size_t maxNodeIndex, size_t timestepIndex = 0) const;

        //! \brief Per boundary-condition-element heat rates from the solved node values.
        //!
        //! Each element's rate is the consistent boundary-side reintegration sum_i((H*u - R)_i)
        //! evaluated with the current node values of the given state variable: the exact heat the
        //! boundary condition exchanges with the discrete solution (the same quantity legacy
        //! Conrad's "eb" routines report per segment). Boundary segments without a boundary
        //! condition (natural/adiabatic) exchange nothing and produce no entry. Nonlinear boundary
        //! conditions evaluate their H/R at the current node state, so on a converged solution the
        //! rates are consistent with the final assembly.
        //! @param variable State variable of the owning domain (Variable::temperature for thermal).
        //! @param timestepIndex Timestep index; zero for steady state.
        [[nodiscard]] std::vector<BoundaryHeatRate> heatRates(Variable variable,
                                                              size_t timestepIndex = 0) const;

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

    private:
        //! Iterates over all active BCs (steady-state + transient at timestepIndex)
        //! and calls func(const IBCLinear2D &) for each.
        template<typename Func>
        void forEachActiveBC(size_t timestepIndex, Func func) const
        {
            for(const auto & bcPtr : m_BCs)
            {
                func(*bcPtr);
            }
            for(const auto & bcTimesteps : m_TransientBCs)
            {
                if(bcTimesteps.size() < timestepIndex)
                {
                    throw std::runtime_error(
                      "Number of boundary conditions provided is less then "
                      "number of timesteps requested.");
                }
                func(*bcTimesteps[timestepIndex]);
            }
        }
    };

}   // namespace HygroThermFEM
