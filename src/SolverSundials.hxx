#pragma once

#include <vector>
#include <memory>

#include "SingleDomain.hxx"
#include "SolutionVariables.hxx"
#include "Solver.hxx"


namespace Sundials
{
    struct SolverPimpl;

    class TransientSingleDomainSundials : public HygroThermFEM::SingleDomainTransientSolver
    {
    public:
        explicit TransientSingleDomainSundials(HygroThermFEM::SingleDomain & domain);

        HygroThermFEM::SingleTimestepSolution
          transient(const std::vector<double> & previousTimestepValues,
                    double t_DTime,
                    size_t timestepIndex) override;

        HygroThermFEM::SingleTimestepSolution
          transient(const std::vector<double> & previousTimestepValues, double t_DTime);

        void setTolerance(double t_Tolerance);

    private:
        std::shared_ptr<SolverPimpl> &
          getSolverPimpl(HygroThermFEM::SingleDomain & domain,
                         const std::vector<double> & previousTimestepValues);

        std::shared_ptr<SolverPimpl> m_pimpl;
    };

    class SolverIDA : public HygroThermFEM::TransientSolver
    {
    public:
        explicit SolverIDA(HygroThermFEM::MultiDomain & domain);

        std::unique_ptr<HygroThermFEM::SingleDomainTransientSolver>
          createSolver(HygroThermFEM::SingleDomain & domain) override;
    };

    // Put vector for now, but this is not going to work in long term because we need to switch
    // between different types of domains.
    std::vector<HygroThermFEM::SingleTimestepSolution>
      transient(HygroThermFEM::SingleDomain & domain,
                const std::vector<double> & previousTimestepValues,
                double t_DTime,
                size_t nTimesteps = 0u);
}   // namespace Sundials