#pragma once

#include <iosfwd>
#include <memory>

#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"
#include "GasCavities.hxx"
#include "BoundaryConditionCoefficients.hxx"
#include "TimestepNotifier.hxx"
#include "TimestepObserver.hxx"
#include "Materials.hxx"
#include "Nodes.hxx"

namespace HygroThermFEM
{
    //! \brief Class to hold solution from single timestep.
    struct SingleSolution
    {
        std::vector<double> solution;   //!< Solution from single transient step
        double dTime;   //!< Timestep for which solution has been performed. Engine can adopt new
                        //!< timestep for which solution will converge.
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
        SingleSolution transient(const std::vector<double> & currentStateValues,
                                 double t_DTime,
                                 size_t timestepIndex = 0);

        //! Returns flux in x and y direction
        [[nodiscard]] std::vector<NodeFlux> flux() const;

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
        void clearModel();

        //! \brief Attach a diagnostic stream for solver iteration logging.
        //! When non-null, the solver prints per-iteration details (NR corrections,
        //! convergence metrics, sub-timestep divisions) to this stream.
        //! Pass nullptr to disable. The caller owns the stream lifetime.
        void setDiagnosticStream(std::ostream * stream);

    protected:
        friend class MultiDomain;

        //! Forms left hand side matrix in steady state solution.
        SquareMatrix steadyStateLeftHandSide();

        //! Form right hand side vector in stead state solution.
        [[nodiscard]] std::vector<double> steadyStateRightHandSide() const;

        //! \brief Forms mass, conductance and H (from boundary condition) matrices.
        //! @param t_DTime Time between two timestep for which calculates are being performed
        //! @param timestepIndex Timestep index used in case of variable timestep input boundary
        //! conditions.
        SquareMatrix transientM_K_H_Matrix(double t_DTime, size_t timestepIndex);

        //! \brief This function retrieves M*U+R vector (where U is state variable)
        //! @param t_PreviousSolution Solution from previous timestep
        //! @param t_DTime Time between two timestep for which calculates are being performed
        //! @param timestepIndex Timestep index used in case of variable timestep input boundary
        //! conditions.
        std::vector<double> transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                                 double t_DTime,
                                                 size_t timestepIndex);

        //! Returns if domain problem is linear.
        [[nodiscard]] bool isLinear() const;

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

        //! \brief Calling timestep calculations
        //! @param currentStateValues Current state values from previous timestep
        //! @param t_DTime Time different for between timesteps
        //! @param timestepIndex Current timestep index used in variable boundary conditions
        std::pair<std::vector<double>, bool> transientTimestep(
          const std::vector<double> & currentStateValues, double t_DTime, size_t timestepIndex);

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
        //! Returns the accepted trial state with reassembled matA/vecB.
        LineSearchResult backtrackingLineSearch(const std::vector<double> & currentSolution,
                                                const std::vector<double> & correctionDU,
                                                double initialRelax,
                                                double currentResidualNorm,
                                                const std::vector<double> & currentStateValues,
                                                double dTime,
                                                size_t timestepIndex);

        //! Updates node values with solution. Implemented by derived classes.
        virtual void updateNodes(const std::vector<double> & solution,
                                 bool updatePreviousTimestep) = 0;

        Nodes & m_NodePool;
        Materials & m_MaterialPool;
        ElementsLinear2D m_Elements;
        BoundaryConditions2D m_BCs;

        FenestrationCommon::GravityVector m_GravityVector{0, -1, 0};

        //! Storage for gas cavities recalculation
        std::unique_ptr<EquivalentGasCavities> gasCavities;

        // Indicates if transient timestep will automatically update previous timestep solution.
        // This should be turned off if used in multidomain because previous timestep should
        // remain constant during iterations.
        bool m_AutomaticUpdatePreviousTimestep;

        //! True when the last successful NR solve had all DOF corrections
        //! clamped at physical bounds (e.g. humidity pinned at 0 or 1).
        //! Indicates re-solving would produce the identical result.
        bool m_LastSolveAtPhysicalBound{false};

        //! Returns whether the last transient solve ended with all DOFs
        //! clamped at physical bounds by limitIncrement.
        [[nodiscard]] bool lastSolveAtPhysicalBound() const;

        //! Optional diagnostic output stream (nullptr = disabled).
        std::ostream * m_DiagStream{nullptr};
    };

}   // namespace HygroThermFEM
