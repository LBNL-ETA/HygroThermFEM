#pragma once

#include <memory>
#include <vector>

#include "SquareMatrix.hxx"

namespace HygroThermFEM
{
    //! Sparse direct solver with symbolic-factorization reuse. A domain's system matrix keeps
    //! the same sparsity pattern for the whole simulation (mesh connectivity is fixed), so the
    //! symbolic analysis runs only when the pattern actually changes -- guarded by matrix size
    //! and non-zero count -- and every other solve pays just the numeric factorization.
    class CLinearSolver
    {
    public:
        CLinearSolver();
        ~CLinearSolver();
        CLinearSolver(CLinearSolver && other) noexcept;
        CLinearSolver & operator=(CLinearSolver && other) noexcept;

        //! Solves lhsMatrix * x = rhsVector, reusing the cached symbolic factorization when the
        //! sparsity pattern matches the previous call.
        std::vector<double> solve(const SquareMatrix & lhsMatrix,
                                  const std::vector<double> & rhsVector);

        //! One-shot solve without cached state (fresh symbolic + numeric factorization).
        static std::vector<double> solveEigen(const SquareMatrix & t_MatrixA,
                                              const std::vector<double> & t_VectorB);

        //! Total symbolic factorizations performed by all solvers since process start.
        //!
        //! Mesh connectivity is fixed for a run, so this is expected to stay in the single
        //! digits. A count that scales with the number of solves means the reuse guard is
        //! missing and every solve is paying for a symbolic analysis it should not need.
        [[nodiscard]] static std::size_t patternAnalysisCount();

        //! Total solves performed by all solvers since process start. Divided by the number of
        //! timesteps this gives the average Newton-Raphson iteration count, which is the figure
        //! that decides whether a run is bound by per-iteration cost or by iteration count.
        [[nodiscard]] static std::size_t solveCount();

    private:
        struct SolverState;
        std::unique_ptr<SolverState> m_State;
    };
}   // namespace HygroThermFEM
