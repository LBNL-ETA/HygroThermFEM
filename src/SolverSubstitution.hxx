#pragma once

#include <vector>

#include "Element2D.hxx"
#include "SolutionVariables.hxx"
#include "Solver.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;
    struct SingleDomain;
}   // namespace HygroThermFEM

namespace HygroThermFEM
{
    //! Calculates steady state for given data
    //! @param domain Domain for which steady state is being calculated
    //! @param output Domain solution
    std::vector<double> steadyState(SingleDomain & domain);

    //! Calculates steady state solution for multiple domains.
    Solution steadyState(HygroThermFEM::MultiDomain & domain);

    class TransientSubstitutionSolver : public TransientSolver
    {
    public:
        // This brings both versions of transient from TransientSolver into scope
        using TransientSolver::transient;

        //! \brief Calculates next timestep values from current values.
        //! @param domain Domain for which transient solution is being calculated
        //! @param previousTimestepValues Current values of state variable or initial condition
        //! @param t_DTime Timestep in transient solution
        //! @param timestepIndex Index for current timestep. Used in variable boundary conditions
        //! case. It is defaulted to zero in case of non-variable boundary condition calculations
        //! are requested.
        //! @return Solution from single transient step. SingleTimestepSolution contains solution
        //! and timestep for which solution has been performed. Engine can adopt new timestep for
        //! which solution will converge.
        SingleTimestepSolution transient(SingleDomain & domain,
                                         const std::vector<double> & previousTimestepValues,
                                         double t_DTime,
                                         size_t timestepIndex) override;

        //! \brief Overriden functions are not allowed to have a default arguments which is the
        //! reason why this function is defined separately.
        SingleTimestepSolution transient(SingleDomain & domain,
                                         const std::vector<double> & previousTimestepValues,
                                         double t_DTime);
    };
}   // namespace HygroThermFEM::Substitution