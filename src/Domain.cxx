#include <cmath>

#include "Domain.hxx"
#include "FEMunique.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"
#include "Nodes.hxx"
#include "TimestepData.hxx"
#include "Materials.hxx"

namespace HygroThermFEM
{
    SquareMatrix IDomain::steadyStateLeftHandSide()
    {
        const auto maxNodeIndex = m_NodePool.maxIndex();
        auto condMat = m_Elements.conductanceMatrix(maxNodeIndex);
        const auto boundaryHMatrix = m_BCs.HMatrix(maxNodeIndex);
        condMat += boundaryHMatrix;

        return condMat;
    }

    std::vector<double> IDomain::steadyStateRightHandSide() const
    {
        return m_BCs.RVector(m_NodePool.maxIndex());
    }

    SquareMatrix IDomain::transientM_K_H_Matrix(const double t_DTime, const size_t timestepIndex)
    {
        const auto maxNodeIndex = m_NodePool.maxIndex();
        const auto MassVec = m_Elements.getLumpedMass(maxNodeIndex, t_DTime);
        auto M_K_H = m_Elements.conductanceMatrix(maxNodeIndex);
        M_K_H = M_K_H.addDiagonal(MassVec);
        M_K_H += m_BCs.HMatrix(maxNodeIndex, timestepIndex);

        return M_K_H;
    }

    std::vector<double>
      IDomain::transientMT_R_Vector(const std::vector<double> & t_PreviousSolution,
                                    const double t_DTime,
                                    const size_t timestepIndex)
    {
        const auto maxNodeIndex = m_NodePool.maxIndex();
        const std::vector<double> MassVec{m_Elements.getLumpedMass(maxNodeIndex, t_DTime)};
        const auto vecR = m_BCs.RVector(maxNodeIndex, timestepIndex)
                          + m_Elements.RVector(maxNodeIndex);

        auto vecB = t_PreviousSolution * MassVec + vecR;

        return vecB;
    }

    std::vector<double> IDomain::steadyState()
    {
        const auto B = steadyStateRightHandSide();
        const auto A = steadyStateLeftHandSide();
        // const auto test = A.toVector();
        return CLinearSolver::solveEigen(A, B);
    }

    SingleSolution IDomain::transient(const std::vector<double> & currentStateValues,
                                      const double t_DTime,
                                      const size_t timestepIndex)
    {
        std::vector<double> solution;
        bool converged{false};
        auto currentDivisionLevel{0u};
        auto maxDivisionLevel{Timesteps::Settings::Instance().getMaxDivisions()};
        double currentDTime{t_DTime};
        double totalTime{0};
        auto stateVariables{currentStateValues};
        unsigned numberOfSubtimesteps{Timesteps::Settings::Instance().getNumberOfSubtimesteps()};

        // In case program failed to converge, it will cut down step to smaller one and will perform
        // multiple consecutive simulations in order to achieve solution at requested timestep.
        while(totalTime < t_DTime)
        {
            notify(currentDivisionLevel, unsigned(totalTime / currentDTime));
            const auto current = transientTimestep(stateVariables, currentDTime, timestepIndex);
            solution = current.first;
            converged = current.second;
            if(!converged)
            {
                currentDTime = currentDTime / numberOfSubtimesteps;
                ++currentDivisionLevel;
                if(currentDivisionLevel > maxDivisionLevel)
                {
                    throw std::runtime_error("Solution failed to converge.");
                }
            }
            else
            {
                stateVariables = solution;
                totalTime += currentDTime;

                // When all DOFs are pinned at physical bounds the NR correction is
                // zero.  Every subsequent sub-timestep would start from the same
                // clamped state, assemble the same system, and produce the same
                // zero correction — so we can skip straight to the end.
                if(m_LastSolveAtPhysicalBound)
                {
                    totalTime = t_DTime;
                }
            }
        }

        return {solution, t_DTime};
    }

    std::pair<std::vector<double>, bool>
      IDomain::transientTimestep(const std::vector<double> & currentStateValues,
                                 const double t_DTime,
                                 const size_t timestepIndex)
    {
        const auto RelaxParameter = SimulationProperties::Instance().relaxationParamter();
        const auto ConvergenceError = SimulationProperties::Instance().errorTolerance();
        const auto MaxIterations = SimulationProperties::Instance().maxNumberOfIterations();

        auto A = transientM_K_H_Matrix(t_DTime, timestepIndex);

        auto B = transientMT_R_Vector(currentStateValues, t_DTime, timestepIndex);

        std::vector<double> solution;
        bool converged{false};
        bool stopIterations{false};

        if(isLinear())
        {
            solution = CLinearSolver::solveEigen(A, B);
            postProcess(solution);
            converged = true;
            m_LastSolveAtPhysicalBound = false;
        }
        else
        {
            // Nonlinear Newton-Raphson iteration with relaxation.
            // Solves the residual equation  A*U = B  iteratively, where A and B
            // depend on the current solution (material properties are state-dependent).
            // Each iteration computes a correction dU = A^{-1} * (B - A*U), applies it
            // with a relaxation parameter, then reassembles A and B at the updated state.
            solution = currentStateValues;
            std::vector<double> normSolution{currentStateValues};

            auto currentNorm = norm(solution);

            size_t numOfIterations = 0;
            double prevMetric = std::numeric_limits<double>::max();
            constexpr size_t minIterationsForOscillationCheck = 4;
            bool convergedViaClamp = false;

            while(!stopIterations && !converged)
            {
                const double previousNorm = currentNorm;
                const auto prevIterSolution = solution;

                // Compute the residual: r = B - A * U_current
                auto temp = A * solution;
                temp = B - temp;

                // Solve for the Newton-Raphson correction: A * dU = r
                auto dU = CLinearSolver::solveEigen(A, temp);

                // Limit the correction per-DOF so the projected solution stays within
                // physical bounds (e.g., humidity in [0, 1]). This prevents overshoot
                // that would cause clamping-induced oscillations in the solver.
                const bool allClamped = limitIncrement(solution, dU, RelaxParameter);

                // Apply the relaxed correction to get the new solution estimate.
                // normSolution uses the full (unrelaxed) correction for convergence check.
                solution = solution + dU * RelaxParameter;
                normSolution = solution + dU;

                postProcess(solution);
                postProcess(normSolution);

                currentNorm = norm(normSolution);

                ++numOfIterations;

                // Update node properties so the next matrix assembly reflects the new state
                updateNodes(solution, m_AutomaticUpdatePreviousTimestep);

                // Reassemble system matrices at the updated state (material properties,
                // boundary conditions, etc. may have changed with the new solution)
                A = transientM_K_H_Matrix(t_DTime, timestepIndex);
                B = transientMT_R_Vector(currentStateValues, t_DTime, timestepIndex);

                // Convergence check: relative change in the solution norm between iterations
                const auto metric = std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12);
                converged = metric <= ConvergenceError;

                // If limitIncrement clamped ALL corrections to (near-)zero, the
                // solution is pinned at physical bounds and cannot improve further.
                if(!converged && allClamped)
                {
                    converged = true;
                    convergedViaClamp = true;
                }

                // Detect Newton-Raphson 2-cycle oscillation.
                // Piecewise-linear material data (e.g., sorption curves with steep kinks)
                // can cause the Jacobian to alternate between two states on consecutive
                // iterations. The convergence metric then stagnates — it stays nearly the
                // same each iteration without decreasing. When the metric changes by less
                // than 1% between consecutive iterations (after a minimum warmup period),
                // the solver is stuck in a 2-cycle. Resolve by averaging the last two
                // iterates, which places the solution at the midpoint of the oscillation
                // and yields a physically reasonable result.
                if(!converged
                   && numOfIterations >= minIterationsForOscillationCheck
                   && std::abs(metric - prevMetric) / (prevMetric + 1e-12) < 0.01)
                {
                    for(size_t idx = 0; idx < solution.size(); ++idx)
                    {
                        solution[idx] = 0.5 * (solution[idx] + prevIterSolution[idx]);
                    }
                    postProcess(solution);
                    converged = true;
                }

                prevMetric = metric;
                stopIterations = numOfIterations > MaxIterations;
            }

            m_LastSolveAtPhysicalBound = convergedViaClamp;
        }

        updateNodes(solution, m_AutomaticUpdatePreviousTimestep);

        return std::make_pair(solution, converged);
    }

    bool IDomain::isLinear() const
    {
        return m_BCs.isLinear() && m_Elements.isLinear();
    }

    IDomain::IDomain(Nodes & nodePool,
                     Materials & materialPool,
                     const bool automaticUpdateOfPreviousTimestep) :
        m_NodePool(nodePool),
        m_MaterialPool(materialPool),
        gasCavities(nullptr),
        m_AutomaticUpdatePreviousTimestep(automaticUpdateOfPreviousTimestep)
    {}

    std::vector<NodeFlux> IDomain::flux() const
    {
        return m_Elements.flux(m_NodePool.maxIndex());
    }

    bool IDomain::lastSolveAtPhysicalBound() const
    {
        return m_LastSolveAtPhysicalBound;
    }

    bool IDomain::limitIncrement(const std::vector<double> & /*currentSolution*/,
                                 std::vector<double> & /*increment*/,
                                 const double /*relaxParameter*/) const
    {
        // Default: no limiting. Derived classes override for physical bounds.
        return false;
    }

    void IDomain::postProcess(std::vector<double> &)
    {
        if(gasCavities == nullptr)
        {
            gasCavities = std::make_unique<EquivalentGasCavities>(m_NodePool, m_MaterialPool, m_Elements);
            gasCavities->setGravityVector(m_GravityVector);
        }
        gasCavities->update();
    }

    void IDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        m_GravityVector = gravityVector;
        if(gasCavities != nullptr)
        {
            gasCavities->setGravityVector(gravityVector);
        }
    }

    void IDomain::clearModel()
    {
        m_NodePool.clear();
        m_MaterialPool.clear();
        m_BCs.clear();
        m_Elements.clearElements();
    }        

}   // namespace HygroThermFEM
