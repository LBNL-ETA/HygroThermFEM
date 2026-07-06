#pragma once

#include <iosfwd>
#include <limits>
#include <memory>
#include <optional>

#include "lbnl/expected.hxx"

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "GasCavities.hxx"
#include "BoundaryConditionCoefficients.hxx"
#include "TimestepNotifier.hxx"
#include "TimestepObserver.hxx"
#include "Materials.hxx"
#include "Nodes.hxx"
#include "SolverSettings.hxx"

namespace HygroThermFEM
{
    //! \brief Class to hold solution from single timestep.
    struct SingleSolution
    {
        std::vector<double> solution;   //!< Solution from single transient step
        double dTime;   //!< Timestep for which solution has been performed. Engine can adopt new
                        //!< timestep for which solution will converge.
    };

    //! \brief Reason a transient solve failed to converge.
    //!
    //! Non-convergence is an expected outcome (not a programmer error), so it is reported as the
    //! error arm of an ExpectedExt rather than thrown from the hot path.
    struct SolverError
    {
        double dtTried;        //!< Timestep at which the final (failed) attempt was made
        size_t divisions;      //!< Number of subdivision levels exhausted before giving up
        double lastResidual;   //!< Residual norm of the last (failed) attempt
    };

    //! \brief Interface that will keep all elements and boundary conditions together.
    //!
    //! One domain will solve single differential equation and therefore, single domain will
    //! represent thermal, moisture or pressure separately.
    class IDomain : public Timesteps::TimestepNotifier
    {
    public:
        virtual ~IDomain() = default;
        //! Domain construction.
        explicit IDomain(
          Nodes & nodePool,           //!< Reference to the NodePool for node lookups
          Materials & materialPool,   //!< Reference to the MaterialPool for material lookups
          bool automaticUpdateOfPreviousTimestep =
            true   //!< When solver finds solution, previous timestep will be automatically updated
                   //!< by default. This should be set to false is Domain is used in outside
                   //!< iterations for solving more complex problems in which case previous timestep
                   //!< should not be updated till convergence is achieved.
        );

        //! Calculates steady state for given data
        std::vector<double> steadyState();

        //! \brief Calculates next timestep values from current (initial) values
        //! @param currentStateValues Current values of state variable or initial condition
        //! @param t_DTime Timestep in transient solution
        //! @param timestepIndex Index for current timestep. Used in variable boundary conditions
        //! case. It is defaulted to zero in case of non-variable boundary condition calculations
        //! are requested.
        //! Returns the solved timestep, or a SolverError if convergence could not be reached even
        //! after exhausting sub-timestep subdivision. Non-convergence is reported, never thrown.
        lbnl::ExpectedExt<SingleSolution, SolverError> transient(
          const std::vector<double> & currentStateValues, double t_DTime, size_t timestepIndex = 0);

        //! Runs a multi-step transient simulation on a single domain,
        //! collecting the solution vector at each timestep.
        //! @param variable Which node property to use as initial condition
        //! @param dTime Timestep duration
        //! @param numSteps Number of timesteps to run
        std::vector<std::vector<double>> transientMultiStep(Variable variable,
                                                            double dTime,
                                                            size_t numSteps);

        //! Returns flux in x and y direction
        [[nodiscard]] std::vector<NodeFlux> flux() const;

        //! \brief Per boundary-condition-element heat rates from the solved node values.
        //! The consistent boundary-side reintegration sum((H*u - R)) per boundary element;
        //! positive when heat leaves the domain (see BoundaryConditions2D::heatRates).
        //! @param variable State variable of this domain (Variable::temperature for thermal).
        [[nodiscard]] std::vector<BoundaryHeatRate> boundaryHeatRates(Variable variable) const;

        //! \brief Read-only access to the domain's elements.
        //! Used by the error-estimator adapter to assemble per-element Gauss-point data from a
        //! solved domain.
        [[nodiscard]] const ElementsLinear2D & elements() const;

        //! Adds element into domain
        virtual void createElement(
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< SolidMaterial that will be assigned to the element
          ) = 0;

        //! \brief Sets new gravity vector and performs new calculations
        //!
        //! @param gravityVector Direction of gravity
        void setGravityVector(const FenestrationCommon::GravityVector & gravityVector);

        //! \brief Deletes geometry and clears up boundary conditions.
        //! Virtual so derived domains can also clear domain-specific state (e.g. ThermalDomain's
        //! enclosure-radiation coordinators) between consecutive runs on the same domain.
        virtual void clearModel();

        //! \brief Attach a diagnostic stream for solver iteration logging.
        //! When non-null, the solver prints per-iteration details (NR corrections,
        //! convergence metrics, sub-timestep divisions) to this stream.
        //! Pass nullptr to disable. The caller owns the stream lifetime.
        void setDiagnosticStream(std::ostream * stream);

        //! \brief Inject solver configuration for this domain.
        //!
        //! When not set, the domain reads the process-global SimulationProperties / Timesteps
        //! singletons at solve time, preserving the historical behaviour.
        void setSolverSettings(const SolverSettings & settings);

    protected:
        //! Some domains require post-processing of results. Good example is
        //! moisture domain where humidity cannot go over 1.0 or lower than one.
        //! With certain set of boundary conditions and long enough time-step,
        //! solution can achieve such state and post processing should prevent it.
        virtual void postProcess(std::vector<double> & solution);

        //! Limits Newton-Raphson increment per DOF so the projected solution
        //! stays within physical bounds. Called before applying dU.
        virtual bool limitIncrement(const std::vector<double> & currentSolution,
                                    std::vector<double> & increment,
                                    double relaxParameter) const;

        //! Updates node values with solution. Implemented by derived classes.
        virtual void updateNodes(const std::vector<double> & solution,
                                 bool updatePreviousTimestep) = 0;

        Nodes & m_NodePool;
        Materials & m_MaterialPool;
        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;

    private:
        //! Forms left hand side matrix in steady state solution.
        SquareMatrix steadyStateLeftHandSide();

        //! Form right hand side vector in steady state solution.
        [[nodiscard]] std::vector<double> steadyStateRightHandSide() const;

        //! \brief Forms mass, conductance and H (from boundary condition) matrices.
        SquareMatrix transientM_K_H_Matrix(double t_DTime, size_t timestepIndex);

        //! \brief This function retrieves M*U+R vector (where U is state variable)
        std::vector<double> transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                                 double t_DTime,
                                                 size_t timestepIndex);

        //! Returns if domain problem is linear.
        [[nodiscard]] bool isLinear() const;

        //! Returns the injected solver settings, or a snapshot of the global singletons when none
        //! has been injected (preserving the historical live-read behaviour).
        [[nodiscard]] SolverSettings solverSettings() const;

        //! Result of a single timestep solve.
        struct TimestepResult
        {
            std::vector<double> solution;
            bool converged;
            double residualNorm = 0.0;   //!< Residual norm of this attempt (used to report failure)
        };

        //! Tracks adaptive damping across Newton-Raphson iterations.
        struct AdaptiveDamping
        {
            double prevDuNorm = (std::numeric_limits<double>::max)();
            double factor = 1.0;
            size_t consecutiveGrowth = 0;

            void update(double rawDuNorm, size_t iteration, size_t minIteration);
        };

        //! Mutable state carried across Newton-Raphson iterations.
        struct NRLoopState
        {
            std::vector<double> solution;
            SquareMatrix matA;
            std::vector<double> vecB;
            double currentNorm;
            double prevMetric = (std::numeric_limits<double>::max)();
            size_t numOfIterations = 0;
            bool convergedViaClamp = false;
            bool converged = false;
            AdaptiveDamping damping;
        };

        //! \brief Calling timestep calculations
        TimestepResult transientTimestep(
          const std::vector<double> & currentStateValues, double t_DTime, size_t timestepIndex);

        //! Performs one Newton-Raphson iteration: correction, clamp check,
        //! line search, convergence and oscillation checks.
        void performNRIteration(NRLoopState & state,
                                const std::vector<double> & currentStateValues,
                                double relaxParameter,
                                double convergenceError,
                                size_t maxIterations,
                                double dTime,
                                size_t timestepIndex);

        //! Result of a backtracking line search step.
        struct LineSearchResult
        {
            std::vector<double> solution;
            SquareMatrix matA;
            std::vector<double> vecB;
            double effectiveRelax;
        };

        //! Backtracking line search (Armijo-style). Shrinks the relaxation
        //! factor until the residual norm decreases or attempts are exhausted.
        LineSearchResult backtrackingLineSearch(const std::vector<double> & currentSolution,
                                                const std::vector<double> & correctionDU,
                                                double initialRelax,
                                                double currentResidualNorm,
                                                const std::vector<double> & currentStateValues,
                                                double dTime,
                                                size_t timestepIndex);

        FenestrationCommon::GravityVector m_GravityVector{0, -1, 0};

        //! Storage for gas cavities recalculation. Empty until the first post-process builds it
        //! lazily; std::optional makes the "already created?" check and ownership explicit without
        //! a heap allocation.
        std::optional<EquivalentGasCavities> gasCavities;

        //! Optional injected solver configuration. Empty means "read the global singletons live".
        std::optional<SolverSettings> m_SolverSettings;

        bool m_AutomaticUpdatePreviousTimestep;
        bool m_LastSolveAtPhysicalBound{false};

        [[nodiscard]] bool lastSolveAtPhysicalBound() const;

        std::ostream * m_DiagStream{nullptr};
    };

}   // namespace HygroThermFEM
