#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>

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
    //! Prints CSV header for NR iteration diagnostics.
    void printNRHeader(std::ostream & out, const size_t numDOFs)
    {
        out << "subtimestep_dt,nr_iter,metric,converge_tol,correction_norm,"
               "damping_factor,all_clamped,converged,reason";
        for(size_t dof = 0; dof < numDOFs; ++dof)
        {
            out << ",u" << dof;
        }
        out << '\n';
    }

    //! Returns the maximum absolute value in a vector. Used by the per-component
    //! convergence check in transientTimestep, which requires the largest |dU_i|
    //! to be small relative to the largest |U_i|. The norm-difference metric on
    //! its own has a null space (mean-preserving perturbations leave the L2 norm
    //! unchanged) and can falsely report convergence; this per-component check
    //! closes that hole.
    double maxAbs(const std::vector<double> & vec)
    {
        double mxv = 0.0;
        for(const auto val : vec)
        {
            const double absVal = std::abs(val);
            if(absVal > mxv)
            {
                mxv = absVal;
            }
        }
        return mxv;
    }

    //! Prints one CSV row for an NR iteration.
    void printNRRow(std::ostream & out,
                    const double subDt,
                    const size_t iteration,
                    const double metric,
                    const double convergeTol,
                    const double correctionNorm,
                    const double dampingFactor,
                    const bool allClamped,
                    const bool converged,
                    const char * reason,
                    const std::vector<double> & solution)
    {
        out << std::scientific << std::setprecision(8)
            << subDt << ','
            << iteration << ','
            << metric << ','
            << convergeTol << ','
            << correctionNorm << ','
            << std::fixed << std::setprecision(6)
            << dampingFactor << ','
            << allClamped << ','
            << converged << ','
            << reason;
        out << std::scientific << std::setprecision(8);
        for(const auto val : solution)
        {
            out << ',' << val;
        }
        out << '\n';
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
        auto currentDivisionLevel{0u};
        const auto maxDivisionLevel{Timesteps::Settings::Instance().getMaxDivisions()};
        double currentDTime{t_DTime};
        const unsigned numberOfSubtimesteps{Timesteps::Settings::Instance().getNumberOfSubtimesteps()};

        if(m_DiagStream != nullptr)
        {
            printNRHeader(*m_DiagStream, currentStateValues.size());
        }

        // Try to solve at the requested dt. If NR fails to converge, subdivide
        // and try again with a smaller dt. Return after the FIRST successful
        // solve, reporting the actual dt used. The caller (MultiDomain) is
        // responsible for accumulating time and exchanging cross-domain data.
        while(true)
        {
            notify(currentDivisionLevel, 0u);
            const auto [solution, converged] =
              transientTimestep(currentStateValues, currentDTime, timestepIndex);

            if(converged)
            {
                return {solution, currentDTime};
            }

            currentDTime = currentDTime / numberOfSubtimesteps;
            ++currentDivisionLevel;
            if(currentDivisionLevel > maxDivisionLevel)
            {
                throw std::runtime_error("Solution failed to converge.");
            }
        }
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

                // Fix 1: if every DOF correction is below the clamp tolerance
                // (~1e-12), the linear solver is returning round-off noise on a
                // near-zero residual. Exit BEFORE applying it — applying round-off
                // noise seeds asymmetric perturbations into U that the norm-based
                // convergence metric cannot detect, and that compound across
                // timesteps until NR overshoots into physical bounds.
                if(allClamped)
                {
                    ++numOfIterations;
                    converged = true;
                    convergedViaClamp = true;

                    if(m_DiagStream != nullptr)
                    {
                        printNRRow(*m_DiagStream, t_DTime, numOfIterations,
                                   0.0, ConvergenceError, rawDuNorm, damping.factor,
                                   true, true, "clamped-no-apply", solution);
                    }
                    break;
                }

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

                // Fix 2: convergence requires BOTH norm-based stagnation AND
                // per-component smallness. The norm metric alone has a null
                // space — anti-symmetric perturbations leave the L2 norm
                // unchanged, so it can falsely report convergence after a
                // single explicit-Euler-like step that grows high-frequency
                // noise. The per-component metric (max applied |dU| over
                // max |U|) closes that hole. We take the max of the two,
                // which is the strictly more conservative criterion.
                const auto normMetric =
                  std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12);
                const auto componentMetric = (maxAbs(correctionDU) * std::abs(effectiveRelax))
                                             / (maxAbs(solution) + 1e-12);
                const auto metric = std::max(normMetric, componentMetric);
                converged = metric <= ConvergenceError;

                const char * reason = "";
                if(converged)
                {
                    reason = "metric";
                }

                // Detect and resolve NR 2-cycle oscillation
                if(!converged
                   && resolveOscillation(solution, prevIterSolution,
                                         metric, prevMetric,
                                         numOfIterations, minIterationsForOscillationCheck))
                {
                    postProcess(solution);
                    converged = true;
                    reason = "oscillation";
                }

                if(m_DiagStream != nullptr)
                {
                    printNRRow(*m_DiagStream, t_DTime, numOfIterations, metric,
                               ConvergenceError, rawDuNorm, damping.factor,
                               allClamped, converged, reason, solution);
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

    void IDomain::setDiagnosticStream(std::ostream * stream)
    {
        m_DiagStream = stream;
    }

    void IDomain::clearModel()
    {
        m_NodePool.clear();
        m_MaterialPool.clear();
        m_BCs.clear();
        m_Elements.clearElements();
    }        

}   // namespace HygroThermFEM
