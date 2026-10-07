#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <ostream>

#include "Domain.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"
#include "Nodes.hxx"
#include "Materials.hxx"

namespace HygroThermFEM
{
    void IDomain::AdaptiveDamping::update(const double rawDuNorm,
                                          const size_t iteration,
                                          const size_t minIteration)
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
}   // namespace HygroThermFEM

namespace
{
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

    // Residual-reduction convergence parameters (D3, moisture). A solve is converged once the
    // free-DOF residual has dropped to reductionFactor * (its value at the initial guess). Being
    // relative to each problem's own initial residual, this is settings- and scale-independent.
    // The absolute floor handles the case where the initial guess already (nearly) solves the
    // system, so the reduction target would otherwise be unreachably tiny.
    constexpr double residualReductionFactor = 1e-6;
    constexpr double residualAbsoluteFloor = 1e-12;

    ConvergenceResult evaluateConvergence(
      std::vector<double> & solution,
      const std::vector<double> & prevIterSolution,
      const std::vector<double> & correctionDU,
      const double effectiveRelax,
      const double previousNorm,
      const double currentNorm,
      const double prevMetric,
      const double convergenceError,
      const bool useResidual,
      const double currentFreeResidual,
      const double initialFreeResidual,
      const size_t numOfIterations,
      const size_t minIterationsForOscillationCheck)
    {
        const auto metric = convergenceMetric(
          previousNorm, currentNorm, correctionDU, effectiveRelax, solution);

        // Thermal keeps the historical change-metric criterion; moisture converges when the
        // free-DOF residual has been reduced by residualReductionFactor from the initial guess.
        const bool primaryConverged =
          useResidual
            ? (currentFreeResidual
               <= residualReductionFactor * initialFreeResidual + residualAbsoluteFloor)
            : (metric <= convergenceError);

        if(primaryConverged)
        {
            return {metric, ConvergeReason::metric};
        }

        // The oscillation-averaging early exit is a crutch for the change-metric: it declares
        // convergence at a stagnated 2-cycle midpoint. Under residual-reduction convergence
        // (moisture) that would accept a not-yet-solved, mass-leaking iterate before the
        // residual is driven down, so it is disabled there -- the loop instead runs to the
        // reduction target or to best-effort acceptance of the lowest-residual iterate.
        if(!useResidual
           && resolveOscillation(solution, prevIterSolution, metric, prevMetric,
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
        auto vecR = m_BCs.RVector(m_NodePool.maxIndex());
        vecR += m_Elements.volumetricSourceVector(m_NodePool.maxIndex());
        return vecR;
    }

    IDomain::TransientSystem
      IDomain::transientSystem(const std::vector<double> & t_PreviousSolution,
                               const double t_DTime,
                               const size_t timestepIndex)
    {
        const auto maxNodeIndex = m_NodePool.maxIndex();
        const auto MassVec = m_Elements.getLumpedMass(maxNodeIndex, t_DTime);

        auto M_K_H = m_Elements.conductanceMatrix(maxNodeIndex);
        M_K_H.addToDiagonal(MassVec);
        M_K_H += m_BCs.HMatrix(maxNodeIndex, timestepIndex);

        auto vecR = m_BCs.RVector(maxNodeIndex, timestepIndex);
        vecR += m_Elements.RVector(maxNodeIndex);
        vecR += m_Elements.volumetricSourceVector(maxNodeIndex);
        vecR += t_PreviousSolution * MassVec;

        return {std::move(M_K_H), std::move(vecR)};
    }

    std::vector<std::vector<double>> IDomain::transientMultiStep(const Variable variable,
                                                                  const double dTime,
                                                                  const size_t numSteps)
    {
        auto values = m_NodePool.properties(variable);
        std::vector<std::vector<double>> results;

        for(size_t step = 0; step < numSteps; ++step)
        {
            const auto stepResult = transient(values, dTime);
            if(!stepResult.has_value())
            {
                throw std::runtime_error("Solution failed to converge.");
            }
            values = stepResult.value().solution;
            results.push_back(values);
        }

        return results;
    }

    std::optional<std::pair<double, double>> IDomain::solutionBounds() const
    {
        return std::nullopt;
    }

    namespace
    {
        enum class PinState : char
        {
            none,
            lower,
            upper
        };

        //! Pins degrees of freedom that violate the bounds; returns true when any new pin was
        //! added.
        bool pinViolators(const std::vector<double> & solution,
                          const std::pair<double, double> & bounds,
                          std::vector<PinState> & pins)
        {
            const auto & [lowerBound, upperBound] = bounds;
            bool changed{false};
            for(std::size_t index = 0u; index < solution.size(); ++index)
            {
                if(pins[index] != PinState::none)
                {
                    continue;
                }
                if(solution[index] > upperBound)
                {
                    pins[index] = PinState::upper;
                    changed = true;
                }
                else if(solution[index] < lowerBound)
                {
                    pins[index] = PinState::lower;
                    changed = true;
                }
            }
            return changed;
        }

        //! Releases pinned degrees of freedom whose own (unpenalized) equation pulls them back
        //! inside the bounds -- without this the first Glaser-style supersaturated solve
        //! over-pins whole regions that belong strictly below saturation. Returns true when any
        //! pin was released.
        bool releaseInactivePins(const SquareMatrix & lhsMatrix,
                                 const std::vector<double> & rhsVector,
                                 const std::vector<double> & solution,
                                 const std::pair<double, double> & bounds,
                                 std::vector<PinState> & pins)
        {
            // Margin against add/release chattering right at the bound.
            constexpr double releaseMargin{1e-9};
            const auto & [lowerBound, upperBound] = bounds;
            const auto residual = rhsVector - lhsMatrix * solution;
            bool changed{false};
            for(std::size_t index = 0u; index < solution.size(); ++index)
            {
                if(pins[index] == PinState::none)
                {
                    continue;
                }
                // Where the DOF's original row equation would place it, its neighbours held
                // fixed: u_i + r_i / A_ii.
                const double freeValue{solution[index]
                                       + residual[index] / lhsMatrix(index, index)};
                const bool releaseUpper{pins[index] == PinState::upper
                                        && freeValue < upperBound - releaseMargin};
                const bool releaseLower{pins[index] == PinState::lower
                                        && freeValue > lowerBound + releaseMargin};
                if(releaseUpper || releaseLower)
                {
                    pins[index] = PinState::none;
                    changed = true;
                }
            }
            return changed;
        }

        //! Original system plus a large diagonal penalty on every pinned degree of freedom --
        //! the penalty dominates every physical conductance in the row, turning it into a
        //! Dirichlet condition at the bound.
        std::pair<SquareMatrix, std::vector<double>>
          penalizedSystem(const SquareMatrix & lhsMatrix,
                          const std::vector<double> & rhsVector,
                          const std::pair<double, double> & bounds,
                          const std::vector<PinState> & pins)
        {
            constexpr double penalty{1e12};
            const auto & [lowerBound, upperBound] = bounds;
            auto matrix{lhsMatrix};
            auto rhs{rhsVector};
            for(std::size_t index = 0u; index < pins.size(); ++index)
            {
                if(pins[index] == PinState::none)
                {
                    continue;
                }
                matrix(index, index) += penalty;
                rhs[index] += penalty * (pins[index] == PinState::upper ? upperBound : lowerBound);
            }
            return {std::move(matrix), std::move(rhs)};
        }
    }   // namespace

    std::vector<double> IDomain::steadyState()
    {
        const auto lhsMatrix = steadyStateLeftHandSide();
        const auto rhsVector = steadyStateRightHandSide();
        auto solution = m_LinearSolver.solve(lhsMatrix, rhsVector);

        const auto bounds{solutionBounds()};
        if(!bounds.has_value())
        {
            return solution;
        }

        // Primal active set WITH release: pin bound violators, re-solve, release pins whose
        // equations pull back inside the bounds, repeat until the set is stable. The transient
        // solver's clamped Newton resolves the same complementarity condition; add-only pinning
        // is not equivalent and freezes the first pass's supersaturated region at the bound.
        // The pass cap guards against chattering on degenerate systems.
        std::vector<PinState> pins(solution.size(), PinState::none);
        const std::size_t maxPasses{solution.size() + 10u};
        for(std::size_t pass = 0u; pass < maxPasses; ++pass)
        {
            const bool anyPinned{pinViolators(solution, bounds.value(), pins)};
            const bool anyReleased{
              releaseInactivePins(lhsMatrix, rhsVector, solution, bounds.value(), pins)};
            if(!anyPinned && !anyReleased)
            {
                break;
            }
            const auto [matrix, rhs] = penalizedSystem(lhsMatrix, rhsVector, bounds.value(), pins);
            solution = m_LinearSolver.solve(matrix, rhs);
        }
        return solution;
    }

    lbnl::ExpectedExt<SingleSolution, SolverError> IDomain::transient(
      const std::vector<double> & currentStateValues,
      const double t_DTime,
      const size_t timestepIndex)
    {
        const auto settings = solverSettings();
        auto currentDivisionLevel{0u};
        const auto maxDivisionLevel{settings.maxDivisions};
        double currentDTime{t_DTime};
        const unsigned numberOfSubtimesteps{settings.numberOfSubtimesteps};

        if(m_DiagStream != nullptr)
        {
            printNRHeader(*m_DiagStream, currentStateValues.size());
        }

        // Try to solve at the requested dt. If NR fails to converge, subdivide
        // and try again with a smaller dt. Return after the FIRST successful
        // solve, reporting the actual dt used. The caller (MultiDomain) is
        // responsible for accumulating time and exchanging cross-domain data.
        bool retriedForBoundHit = false;
        while(true)
        {
            notify(currentDivisionLevel, 0u);
            auto result = transientTimestep(currentStateValues, currentDTime, timestepIndex);

            if(result.converged && retryStepOnNewBoundHit() && !retriedForBoundHit
               && hasNewBoundHits(currentStateValues, result.solution))
            {
                // Initial-condition/boundary shock: the accepted solution drove DOFs onto a
                // physical bound they did not start at (overshoot clamped by postProcess).
                // Resolve the shock in time by redoing this one step as a bounded substep
                // chain. Single-shot: if the resolved step still ends at a bound, that is
                // accepted as physical.
                retriedForBoundHit = true;
                result = resolveShockStep(currentStateValues, currentDTime, timestepIndex);
            }

            if(result.converged)
            {
                return SingleSolution{result.solution, currentDTime};
            }

            currentDTime = currentDTime / numberOfSubtimesteps;
            ++currentDivisionLevel;
            if(currentDivisionLevel > maxDivisionLevel)
            {
                return lbnl::Unexpected(
                  SolverError{currentDTime, currentDivisionLevel, result.residualNorm});
            }
        }
    }

    bool IDomain::hasNewBoundHits(const std::vector<double> & startValues,
                                  const std::vector<double> & solution) const
    {
        const auto startConstrained = constrainedDofs(startValues);
        const auto solutionConstrained = constrainedDofs(solution);
        for(std::size_t idx = 0; idx < solution.size(); ++idx)
        {
            if(solutionConstrained[idx] && !startConstrained[idx])
            {
                return true;
            }
        }
        return false;
    }

    IDomain::TimestepResult IDomain::resolveShockStep(const std::vector<double> & startValues,
                                                      const double t_DTime,
                                                      const size_t timestepIndex)
    {
        constexpr unsigned shockSubsteps = 8;
        const auto entryPreviousValues =
          m_NodePool.properties(stateVariable(), Timestep::Previous);
        const double subDTime = t_DTime / shockSubsteps;

        auto values = startValues;
        TimestepResult result{.solution = startValues, .converged = true};
        for(unsigned substep = 0; substep < shockSubsteps; ++substep)
        {
            // Make the nodes' previous-timestep value the substep start, so the secant
            // storage capacity telescopes (conserves) per substep. Double update: the
            // first call sets the current value, the second rolls it into previous.
            updateNodes(values, false);
            updateNodes(values, true);
            result = transientTimestep(values, subDTime, timestepIndex);
            if(!result.converged)
            {
                break;
            }
            values = result.solution;
        }

        // Restore the previous-timestep values owned by the caller (coupling repeats rely
        // on previous == step start), then re-impose the final solution as current.
        updateNodes(entryPreviousValues, false);
        updateNodes(entryPreviousValues, true);
        updateNodes(result.solution, m_AutomaticUpdatePreviousTimestep);
        return result;
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
            auto [matA, vecB] = transientSystem(currentStateValues, dTime, timestepIndex);

            const double trialResidualNorm = norm(vecB - matA * trialSolution);

            if(trialResidualNorm < currentResidualNorm || attempt == maxAttempts - 1)
            {
                return {std::move(trialSolution), std::move(matA),
                        std::move(vecB), effectiveRelax};
            }
            effectiveRelax *= shrinkFactor;
        }

#if defined(_MSC_VER)
        __assume(false);
#elif defined(__GNUC__) || defined(__clang__)
        __builtin_unreachable();
#endif
    }

    double IDomain::freeDofResidualNorm(const std::vector<double> & solution,
                                        const SquareMatrix & matA,
                                        const std::vector<double> & vecB) const
    {
        const auto residual = vecB - matA * solution;
        const auto constrained = constrainedDofs(solution);
        double sumSq = 0.0;
        for(std::size_t i = 0; i < residual.size(); ++i)
        {
            if(!constrained[i])
            {
                sumSq += residual[i] * residual[i];
            }
        }
        return std::sqrt(sumSq);
    }

    void IDomain::performNRIteration(NRLoopState & state,
                                      const std::vector<double> & currentStateValues,
                                      const double relaxParameter,
                                      const double convergenceError,
                                      const size_t maxIterations,
                                      const double dTime,
                                      const size_t timestepIndex)
    {
        constexpr size_t minIterationsForOscillationCheck = 4;
        const double previousNorm = state.currentNorm;
        const auto prevIterSolution = state.solution;

        const auto residual = state.vecB - state.matA * state.solution;
        auto correctionDU = m_LinearSolver.solve(state.matA, residual);
        const double rawDuNorm = norm(correctionDU);

        if(limitIncrement(state.solution, correctionDU, relaxParameter))
        {
            ++state.numOfIterations;
            state.convergedViaClamp = true;
            state.converged = true;
            return;
        }

        // Thermal uses the historical dU-growth damping. Moisture instead adapts damping from
        // the residual trend (below), so its damping factor -- already updated at the end of the
        // previous iteration -- is left as-is here.
        if(!useResidualConvergence())
        {
            state.damping.update(rawDuNorm, state.numOfIterations, maxIterations / 2);
        }

        // Moisture ignores the user relaxation parameter and lets the adaptive damping fully
        // control under-relaxation, so the solve is independent of the relaxation setting;
        // thermal keeps the historical relaxParameter * damping behaviour.
        const double baseRelax = useResidualConvergence() ? 1.0 : relaxParameter;
        auto lsResult = backtrackingLineSearch(
          state.solution, correctionDU,
          baseRelax * state.damping.factor, norm(residual),
          currentStateValues, dTime, timestepIndex);
        state.solution = std::move(lsResult.solution);
        state.matA = std::move(lsResult.matA);
        state.vecB = std::move(lsResult.vecB);

        auto normSolution = state.solution + correctionDU;
        postProcess(normSolution);
        state.currentNorm = norm(normSolution);
        ++state.numOfIterations;

        // Free-DOF residual of the current iterate (matA/vecB were reassembled for it in the
        // line search). Moisture converges when this drops to a fraction of the initial-guess
        // residual (residual-reduction, D3); thermal ignores it.
        const bool useResidual = useResidualConvergence();
        const double currentFreeResidual =
          useResidual ? freeDofResidualNorm(state.solution, state.matA, state.vecB) : 0.0;

        // Track the lowest-residual iterate so best-effort acceptance returns the best point
        // visited (not an oscillating last iterate) when the reduction target isn't met.
        if(useResidual && currentFreeResidual < state.bestResidual)
        {
            state.bestResidual = currentFreeResidual;
            state.bestSolution = state.solution;
        }

        // Adaptive under-relaxation (moisture): if this step increased the free residual the
        // iteration is diverging (typical at relaxation 1.0 when wetting toward saturation, where
        // the sorption capacity is very stiff), so halve the damping factor; if it decreased,
        // relax the damping back toward 1. The reduced factor feeds the next iteration's line
        // search, converting a divergent fixed-relaxation solve into a stable under-relaxed one.
        // Start at full relaxation (damping factor 1.0, fresh each timestep) and ratchet it down
        // when a step increases the free residual (divergence); recover it toward 1.0 when a step
        // makes progress, so a case that only needs damping through an early stiff transient can
        // speed back up once stable. Purely monotonic decrease converges too slowly for cases
        // that recover.
        if(useResidual)
        {
            constexpr double minDamping = 1.0 / 64.0;
            if(currentFreeResidual > state.prevFreeResidual)
            {
                state.damping.factor = std::max(state.damping.factor * 0.5, minDamping);
            }
            else
            {
                state.damping.factor = std::min(state.damping.factor * 1.5, 1.0);
            }
            state.prevFreeResidual = currentFreeResidual;
        }

        const auto result = evaluateConvergence(
          state.solution, prevIterSolution, correctionDU, lsResult.effectiveRelax,
          previousNorm, state.currentNorm, state.prevMetric, convergenceError,
          useResidual, currentFreeResidual, state.initialFreeResidual,
          state.numOfIterations, minIterationsForOscillationCheck);

        if(result.reason == ConvergeReason::oscillation)
        {
            postProcess(state.solution);
        }
        state.converged = result.converged();
        state.prevMetric = result.metric;

        if(m_DiagStream != nullptr)
        {
            printNRRow(*m_DiagStream, dTime, state.numOfIterations, result.metric,
                       convergenceError, rawDuNorm, state.damping.factor,
                       false, state.converged, reasonLabel(result.reason), state.solution);
        }
    }

    IDomain::TimestepResult
      IDomain::transientTimestep(const std::vector<double> & currentStateValues,
                                 const double t_DTime,
                                 const size_t timestepIndex)
    {
        const auto settings = solverSettings();
        const auto relaxParameter = settings.relaxationParameter;
        const auto convergenceError = settings.errorTolerance;
        auto maxIterations = settings.maxNumberOfIterations;

        // With the latent heat of fusion active the thermal domain converges on residual
        // reduction (see ThermalDomain::useResidualConvergence): the secant capacity kink
        // contracts the fixed point slowly, so the iteration needs room to drive the
        // residual down before best-effort acceptance -- the user's usual budget (~25)
        // was tuned for the fast change-metric exits.
        if(useResidualConvergence() && !physicsOptions().excludeLatentHeatOfFusion)
        {
            maxIterations = (std::max)(maxIterations, static_cast<std::size_t>(200));
        }

        auto [matA, vecB] = transientSystem(currentStateValues, t_DTime, timestepIndex);

        if(isLinear())
        {
            auto solution = m_LinearSolver.solve(matA, vecB);
            postProcess(solution);
            m_LastSolveAtPhysicalBound = false;
            updateNodes(solution, m_AutomaticUpdatePreviousTimestep);
            return {solution, true, 0.0};
        }

        // Residual of the initial guess over the free DOFs -- the reference against which the
        // residual-reduction convergence test measures progress (D3, moisture only).
        const double initialFreeResidual =
          useResidualConvergence() ? freeDofResidualNorm(currentStateValues, matA, vecB) : 0.0;

        NRLoopState state{
          .solution = currentStateValues,
          .matA = std::move(matA),
          .vecB = std::move(vecB),
          .currentNorm = norm(currentStateValues),
          .initialFreeResidual = initialFreeResidual,
          .damping = {}};

        while(state.numOfIterations <= maxIterations && !state.converged)
        {
            performNRIteration(state, currentStateValues, relaxParameter,
                               convergenceError, maxIterations, t_DTime, timestepIndex);
        }

        // Best-effort acceptance for residual-reduction domains (moisture): if the iteration
        // budget is exhausted without meeting the reduction target, return the best (lowest
        // free-residual) iterate seen and accept it, rather than reporting non-convergence
        // (which would drive subdivision to a throw) or returning an oscillating last iterate.
        // This mirrors the reference solver, which runs its iterations and returns; the
        // reduction criterion still governs early exit and accuracy for the converging cases.
        if(useResidualConvergence())
        {
            if(!state.converged && !state.bestSolution.empty())
            {
                state.solution = state.bestSolution;
            }
            state.converged = true;
        }

        m_LastSolveAtPhysicalBound = state.convergedViaClamp;
        updateNodes(state.solution, m_AutomaticUpdatePreviousTimestep);

        const auto residualNorm = norm(state.vecB - state.matA * state.solution);
        return {state.solution, state.converged, residualNorm};
    }

    bool IDomain::isLinear() const
    {
        return m_BCs.isLinear() && m_Elements.isLinear();
    }

    void IDomain::setSolverSettings(const SolverSettings & settings)
    {
        m_SolverSettings = settings;
    }

    SolverSettings IDomain::solverSettings() const
    {
        return m_SolverSettings ? *m_SolverSettings : SolverSettings::fromGlobals();
    }

    void IDomain::setPhysicsOptions(const PhysicsOptions & options)
    {
        m_PhysicsOptions = options;
    }

    PhysicsOptions IDomain::physicsOptions() const
    {
        return m_PhysicsOptions ? *m_PhysicsOptions : PhysicsOptions::fromGlobals();
    }

    IDomain::IDomain(Nodes & nodePool,
                     Materials & materialPool,
                     const bool automaticUpdateOfPreviousTimestep) :
        m_NodePool(nodePool),
        m_MaterialPool(materialPool),
        m_AutomaticUpdatePreviousTimestep(automaticUpdateOfPreviousTimestep)
    {}

    std::vector<NodeFlux> IDomain::flux() const
    {
        return m_Elements.flux(m_NodePool.maxIndex());
    }

    std::vector<BoundaryHeatRate> IDomain::boundaryHeatRates(const Variable variable) const
    {
        return m_BCs.heatRates(variable);
    }

    const ElementsLinear2D & IDomain::elements() const
    {
        return m_Elements;
    }

    void IDomain::setVolumetricSource(const double value)
    {
        m_Elements.setVolumetricSource(value);
    }

    void IDomain::setVolumetricSource(const std::vector<double> & perElement)
    {
        m_Elements.setVolumetricSource(perElement);
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
        if(!gasCavities)
        {
            gasCavities.emplace(m_NodePool, m_MaterialPool, m_Elements);
            gasCavities->setGravityVector(m_GravityVector);
        }
        gasCavities->update();
    }

    void IDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        m_GravityVector = gravityVector;
        if(gasCavities)
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
