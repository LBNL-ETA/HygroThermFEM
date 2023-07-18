#pragma once

#include "SolutionVariables.hxx"
#include "FEMMath.hxx"

namespace HygroThermFEM
{
    struct MultiDomain;
    struct SingleDomain;

    struct IterationResult
    {
        double error{std::numeric_limits<double>::max()};
        SingleTimestepSolution solution;
    };

    class SingleDomainTransientSolver
    {
    public:
        explicit SingleDomainTransientSolver(SingleDomain & domain);
        ~SingleDomainTransientSolver() = default;

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
        virtual HygroThermFEM::SingleTimestepSolution
          transient(const std::vector<double> & previousTimestepValues,
                    double t_DTime,
                    size_t timestepIndex) = 0;

        IterationResult performDomainIteration(std::vector<double> & currentVariable,
                                               const std::vector<double> & previousTimestepVariable,
                                               double dTime,
                                               size_t timestepIndex);

    protected:
        SingleDomain & m_Domain;
    };

    class TransientSolver
    {
    public:
        explicit TransientSolver(MultiDomain & domain);

        virtual ~TransientSolver() = default;

        //! \brief Calculates next timestep value from current values
        //! \param previousTimestepTemperature vector of nodal temperatures from previous timestep
        //! \param previousTimestepHumidity vector of nodal humidity from previous timestep
        //! \param t_DTime time between two timestep
        //! \param timestepIndex current timestep index used in variable boundary conditions case
        //! \return Solution from single transient step. SingleTimestepSolution contains solution
        //! and timestep for which solution has been performed. Engine can adopt new timestep for
        //! which solution will converge.
        Solution transient(const std::vector<double> & previousTimestepTemperature,
                           const std::vector<double> & previousTimestepHumidity,
                           double t_DTime,
                           size_t timestepIndex);

        Solution transient(const std::vector<double> & previousTimestepTemperature,
                           const std::vector<double> & previousTimestepHumidity,
                           double t_DTime);

        virtual std::unique_ptr<SingleDomainTransientSolver>
          createSolver(HygroThermFEM::SingleDomain & domain) = 0;

    protected:
        std::unique_ptr<SingleDomainTransientSolver> m_TemperatureSolver{nullptr};
        std::unique_ptr<SingleDomainTransientSolver> m_HumiditySolver{nullptr};

        MultiDomain & m_Domain;
    };

}   // namespace HygroThermFEM