#include "Solver.hxx"

#include "NodePool.hxx"
#include "SingleDomain.hxx"

namespace HygroThermFEM
{

    IterationResult
      TransientSolver::performDomainIteration(SingleDomain & domain,
                                              std::vector<double> & currentVariable,
                                              const std::vector<double> & previousVariable,
                                              double dTime,
                                              size_t timestepIndex)
    {

        {
            auto newValueSolution = transient(domain, previousVariable, dTime, timestepIndex);
            auto newValueError =
              HygroThermFEM::errorNorm(newValueSolution.solution, currentVariable);
            updateNodeValues(newValueSolution.solution, baseVariableOf(domain), false);
            currentVariable = newValueSolution.solution;
            return {newValueError, newValueSolution};
        }
    }
}