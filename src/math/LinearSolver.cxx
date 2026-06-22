#pragma warning(push, 0)
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#pragma warning(pop)

#include <stdexcept>

#include "LinearSolver.hxx"

namespace HygroThermFEM
{
    std::vector<double> CLinearSolver::solveEigen(const SquareMatrix & t_MatrixA,
                                                  const std::vector<double> & t_VectorB)
    {
        using Matrix = Eigen::SparseMatrix<double>;
        using Vector = Eigen::VectorXd;

        const Vector B = Vector::Map(t_VectorB.data(), t_VectorB.size());

        // Factorize once: getSparseMatrix() now returns a const reference, so the sparse
        // matrix is no longer copied per analyzePattern/factorize call.
        const Matrix & matrix = t_MatrixA.getSparseMatrix();
        Eigen::SparseLU<Matrix> solver;
        solver.analyzePattern(matrix);
        solver.factorize(matrix);

        if(solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Linear solver failed to factorize the system matrix.");
        }

        const Vector y = solver.solve(B);

        if(solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Linear solver failed to solve the system.");
        }

        return std::vector<double>(y.data(), y.data() + y.size());
    }

}   // namespace HygroThermFEM
