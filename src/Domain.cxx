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

    //! Computes the combined convergence metric: the stricter of the
    //! norm-based stagnation metric and the per-component metric.
    double convergenceMetric(const double previousNorm,
                             const double currentNorm,
                             const std::vector<double> & correctionDU,
                             const double effectiveRelax,
                             const std::vector<double> & solution)
    {
        const auto normMetric =
          std::abs(previousNorm - currentNorm) / (currentNorm + 1e-12);
        const auto componentMetric =
          (maxAbs(correctionDU) * std::abs(effectiveRelax)) / (maxAbs(solution) + 1e-12);
        return std::max(normMetric, componentMetric);
    }

    enum class ConvergeReason { none, metric, oscillation };

    const char * reasonLabel(const ConvergeReason reason)
    {
        switch(reason)
        {
            case ConvergeReason::metric: return "metric";
            case ConvergeReason::oscillation: return "oscillation";
            default: return "";
        }
    }

    struct ConvergenceResult
    {
        double metric;
        ConvergeReason reason = ConvergeReason::none;

        [[nodiscard]] bool converged() const
        {
            return reason != ConvergeReason::none;
        }
    };

    ConvergenceResult evaluateConvergence(
      std::vector<double> & solution,
      const std::vector<double> & prevIterSolution,
      const std::vector<double> & correctionDU,
      const double effectiveRelax,
      const double previousNorm,
      const double currentNorm,
      const double prevMetric,
      const double convergenceError,
      const size_t numOfIterations,
      const size_t minIterationsForOscillationCheck)
    {
        const auto metric = convergenceMetric(
          previousNorm, currentNorm, correctionDU, effectiveRelax, solution);

        if(metric <= convergenceError)
        {
            return {metric, ConvergeReason::metric};
        }

        if(resolveOscillation(solution, prevIterSolution,
                              metric, prevMetric,
                              numOfIterations, minIterationsForOscillationCheck))
        {
            return {metric, ConvergeReason::oscillation};
        }

        return {metric};
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

    IDomain::LineSearchResult IDomain::backtrackingLineSearch(
      const std::vector<double> & currentSolution,
      const std::vector<double> & correctionDU,
      const double initialRelax,
      const double currentResidualNorm,
      const std::vector<double> & currentStateValues,
      const double dTime,
      const size_t timestepIndex)
    {
        constexpr int maxAttempts = 8;
        constexpr double shrinkFactor = 0.5;
        double effectiveRelax = initialRelax;

        for(int attempt = 0; attempt < maxAttempts; ++attempt)
        {
            auto trialSolution = currentSolution + correctionDU * effectiveRelax;
            postProcess(trialSolution);

            updateNodes(trialSolution, m_AutomaticUpdatePreviousTimestep);
            auto matA = transientM_K_H_Matrix(dTime, timestepIndex);
            auto vecB = transientMT_R_Vector(currentStateValues, dTime, timestepIndex);

            const double trialResidualNorm = norm(vecB - matA * trialSolution);

            if(trialResidualNorm < currentResidualNorm || attempt == maxAttempts - 1)
            {
                return {std::move(trialSolution), std::move(matA),
                        std::move(vecB), effectiveRelax};
            }
            effectiveRelax *= shrinkFactor;
        }

        __assume(false);
    }

    std::pair<std::vector<double>, bool>
      IDomain::transientTimestep(const std::vector<double> & currentStateValues,
                                 const double t_DTime,
                                 const size_t timestepIndex)
    {
        const auto relaxParameter = SimulationProperties::Instance().relaxationParamter();
        const auto convergenceError = SimulationProperties::Instance().errorTolerance();
        const auto maxIterations = SimulationProperties::Instance().maxNumberOfIterations();

        auto matA = transientM_K_H_Matrix(t_DTime, timestepIndex);
        auto vecB = transientMT_R_Vector(currentStateValues, t_DTime, timestepIndex);

        if(isLinear())
        {
            auto solution = CLinearSolver::solveEigen(matA, vecB);
            postProcess(solution);
            m_LastSolveAtPhysicalBound = false;
            updateNodes(solution, m_AutomaticUpdatePreviousTimestep);
            return {solution, true};
        }

        auto solution = currentStateValues;
        auto currentNorm = norm(solution);

        size_t numOfIterations = 0;
        double prevMetric = std::numeric_limits<double>::max();
        constexpr size_t minIterationsForOscillationCheck = 4;
        bool convergedViaClamp = false;
        bool converged = false;
        AdaptiveDamping damping;

        while(numOfIterations <= maxIterations && !converged)
        {
            const double previousNorm = currentNorm;
            const auto prevIterSolution = solution;

            const auto residual = vecB - matA * solution;
            auto correctionDU = CLinearSolver::solveEigen(matA, residual);
            const double rawDuNorm = norm(correctionDU);

            if(limitIncrement(solution, correctionDU, relaxParameter))
            {
                ++numOfIterations;
                convergedViaClamp = true;
                converged = true;
                break;
            }

            damping.update(rawDuNorm, numOfIterations, maxIterations / 2);

            auto lsResult = backtrackingLineSearch(
              solution, correctionDU,
              relaxParameter * damping.factor, norm(residual),
              currentStateValues, t_DTime, timestepIndex);
            solution = std::move(lsResult.solution);
            matA = std::move(lsResult.matA);
            vecB = std::move(lsResult.vecB);
            const double effectiveRelax = lsResult.effectiveRelax;

            auto normSolution = solution + correctionDU;
            postProcess(normSolution);
            currentNorm = norm(normSolution);
            ++numOfIterations;

            const auto result = evaluateConvergence(
              solution, prevIterSolution, correctionDU, effectiveRelax,
              previousNorm, currentNorm, prevMetric, convergenceError,
              numOfIterations, minIterationsForOscillationCheck);

            if(result.reason == ConvergeReason::oscillation)
            {
                postProcess(solution);
            }
            converged = result.converged();

            if(m_DiagStream != nullptr)
            {
                printNRRow(*m_DiagStream, t_DTime, numOfIterations, result.metric,
                           convergenceError, rawDuNorm, damping.factor,
                           false, converged, reasonLabel(result.reason), solution);
            }

            prevMetric = result.metric;
        }

        m_LastSolveAtPhysicalBound = convergedViaClamp;
        updateNodes(solution, m_AutomaticUpdatePreviousTimestep);

        return {solution, converged};
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
