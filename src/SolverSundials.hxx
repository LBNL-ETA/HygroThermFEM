#pragma once

#include <vector>
#include <memory>

#include "SolutionVariables.hxx"
#include "Solver.hxx"

namespace HygroThermFEM
{
    struct SingleDomain;
}

namespace Sundials
{
    struct SolverPimpl;

    class SolverIDA : public HygroThermFEM::TransientSolver
    {
    public:
        ~SolverIDA() override;

        HygroThermFEM::SingleTimestepSolution
          transient(HygroThermFEM::SingleDomain & domain,
                    const std::vector<double> & previousTimestepValues,
                    double t_DTime,
                    size_t timestepIndex) override;

        HygroThermFEM::SingleTimestepSolution
          transient(const std::shared_ptr<SolverPimpl> & pimpl,
                    HygroThermFEM::SingleDomain & domain,
                    const std::vector<double> & previousTimestepValues,
                    double t_DTime,
                    size_t timestepIndex);

        HygroThermFEM::SingleTimestepSolution
          transient(HygroThermFEM::SingleDomain & domain,
                    const std::vector<double> & previousTimestepValues,
                    double t_DTime);

        HygroThermFEM::Solution transient(HygroThermFEM::MultiDomain & domain,
                                          const std::vector<double> & previousTimestepTemperature,
                                          const std::vector<double> & previousTimestepHumidity,
                                          double t_DTime,
                                          size_t timestepIndex) override;

        HygroThermFEM::Solution transient(HygroThermFEM::MultiDomain & domain,
                                          const std::vector<double> & previousTimestepTemperature,
                                          const std::vector<double> & previousTimestepHumidity,
                                          double t_DTime);

    private:
        std::shared_ptr<SolverPimpl> m_pimpl;
    };

    // Put vector for now, but this is not going to work in long term because we need to switch
    // between different types of domains.
    std::vector<HygroThermFEM::SingleTimestepSolution>
      transient(HygroThermFEM::SingleDomain & domain,
                const std::vector<double> & previousTimestepValues,
                double t_DTime,
                size_t nTimesteps = 0u);
}   // namespace Sundials