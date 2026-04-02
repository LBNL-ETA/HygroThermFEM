#include <algorithm>
#include <cmath>
#include <limits>

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

namespace
{
    //! Tracks adaptive damping state across Newton-Raphson iterations.
    //! When the raw correction norm grows explosively on consecutive iterations,
    //! the damping factor is halved to stabilise the solver.
    struct AdaptiveDamping
    {
        double prevDuNorm = std::numeric_limits<double>::max();
        double factor = 1.0;
        size_t consecutiveGrowth = 0;

        void update(const double rawDuNorm, const size_t iteration, const size_t minIteration)
        {
            constexpr double growthThreshold = 1000.0;
            constexpr double minDamping = 1.0 / 16.0;
            constexpr size_t requiredConsecutive = 2;

            if(iteration >= minIteration && rawDuNorm > growthThreshold * prevDuNorm)
            {
                ++consecutiveGrowth;
                if(consecutiveGrowth >= requiredConsecutive)
                {
                    factor = std::max(factor * 0.5, minDamping);
                }
            }
            else
            {
                consecutiveGrowth = 0;
                if(factor < 1.0)
                {
                    factor = std::min(factor * 2.0, 1.0);
                }
            }
            prevDuNorm = rawDuNorm;
        }
    };

    //! Detects NR 2-cycle oscillation by checking if the convergence metric
    //! stagnates (changes by less than 1% between consecutive iterations).
    //! If detected, averages the last two iterates to break the cycle.
    //! Returns true if oscillation was resolved.
    bool resolveOscillation(std::vector<double> & solution,
                            const std::vector<double> & prevIterSolution,
                            const double metric,
                            const double prevMetric,
                            const size_t iteration,
                            const size_t minIterations)
    {
        constexpr double stagnationThreshold = 0.01;
        if(iteration < minIterations)
        {
            return false;
        }
        if(std::abs(metric - prevMetric) / (prevMetric + 1e-12) >= stagnationThreshold)
        {
            return false;
        }

        std::ranges::transform(solution, prevIterSolution, solution.begin(),
            [](const double cur, const double prev) { return 0.5 * (cur + prev); });
        return true;
    }
}   // anonymous namespace

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
            // Nonlinear Newton-Raphson iteration with relaxation and adaptive
            // damping.  Solves the residual equation  A*U = B  iteratively, where
            // A and B depend on the current solution (material properties are
            // state-dependent).  Each iteration computes a correction
            // dU = A^{-1} * (B - A*U), applies it with a relaxation parameter, then
            // reassembles A and B at the updated state.
            solution = currentStateValues;
            std::vector<double> normSolution{currentStateValues};
            auto currentNorm = norm(solution);

            size_t numOfIterations = 0;
            double prevMetric = std::numeric_limits<double>::max();
            constexpr size_t minIterationsForOscillationCheck = 4;
            bool convergedViaClamp = false;
            AdaptiveDamping damping;

            while(!stopIterations && !converged)
            {
                const double previousNorm = currentNorm;
                const auto prevIterSolution = solution;

                // Solve for the Newton-Raphson correction: A * dU = B - A*U
                auto residual = B - A * solution;
                auto correctionDU = CLinearSolver::solveEigen(A, residual);
                const double rawDuNorm = norm(correctionDU);

                // Limit per-DOF so the projected solution stays within physical bounds
                const bool allClamped = limitIncrement(solution, correctionDU, RelaxParameter);

                // Adaptive damping: scale down when correction norm explodes
                damping.update(rawDuNorm, numOfIterations, MaxIterations / 2);
                const double effectiveRelax = RelaxParameter * damping.factor;

                // Apply the relaxed correction; normSolution uses full correction
                // for the convergence check
                solution = solution + correctionDU * effectiveRelax;
                normSolution = solution + correctionDU;

                postProcess(solution);
                postProcess(normSolution);
                currentNorm = norm(normSolution);
                ++numOfIterations;

                // Reassemble system at the updated state
                updateNodes(solution, m_AutomaticUpdatePreviousTimestep);
                A = transientM_K_H_Matrix(t_DTime, timestepIndex);
                B = transientMT_R_Vector(currentStateValues, t_DTime, timestepIndex);

                // Convergence check: relative change in solution norm
                const auto metric = std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12);
                converged = metric <= ConvergenceError;

                // All DOFs pinned at physical bounds — cannot improve further
                if(!converged && allClamped)
                {
                    converged = true;
                    convergedViaClamp = true;
                }

                // Detect and resolve NR 2-cycle oscillation
                if(!converged
                   && resolveOscillation(solution, prevIterSolution,
                                         metric, prevMetric,
                                         numOfIterations, minIterationsForOscillationCheck))
                {
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
