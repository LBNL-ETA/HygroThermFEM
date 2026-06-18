#pragma warning(push, 0)
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#pragma warning(pop)

#include "LinearSolver.hxx"

namespace HygroThermFEM
{
    std::vector<double> CLinearSolver::solveEigen(const SquareMatrix & t_MatrixA,
                                                  const std::vector<double> & t_VectorB)
    {
        using Matrix = Eigen::SparseMatrix<double>;
        using Vector = Eigen::VectorXd;

        const auto size = t_MatrixA.size();

        Vector B(size);

        for(auto j = 0u; j < t_VectorB.size(); ++j)
        {
            B[j] = t_VectorB[j];
        }

        Eigen::SparseLU<Matrix> solver;
        solver.analyzePattern(t_MatrixA.getSparseMatrix());
        solver.factorize(t_MatrixA.getSparseMatrix());
        Vector y = solver.solve(B);

        std::vector<double> solution(y.size());
        for(auto i = 0u; i < size; ++i)
        {
            solution[i] = y[i];
        }

        return solution;
    }

}   // namespace HygroThermFEM
