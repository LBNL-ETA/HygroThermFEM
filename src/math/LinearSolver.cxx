#pragma warning(push, 0)
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#pragma warning(pop)

#include <stdexcept>

#include "LinearSolver.hxx"

namespace HygroThermFEM
{
    struct CLinearSolver::SolverState
    {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        Eigen::Index analyzedSize{-1};
        Eigen::Index analyzedNonZeros{-1};
    };

    CLinearSolver::CLinearSolver() : m_State(std::make_unique<SolverState>())
    {}

    CLinearSolver::~CLinearSolver() = default;

    CLinearSolver::CLinearSolver(CLinearSolver && other) noexcept = default;

    CLinearSolver & CLinearSolver::operator=(CLinearSolver && other) noexcept = default;

    std::vector<double> CLinearSolver::solve(const SquareMatrix & lhsMatrix,
                                             const std::vector<double> & rhsVector)
    {
        const Eigen::SparseMatrix<double> & matrix = lhsMatrix.getSparseMatrix();

        const bool samePattern = matrix.rows() == m_State->analyzedSize
                                 && matrix.nonZeros() == m_State->analyzedNonZeros;
        if(!samePattern)
        {
            m_State->solver.analyzePattern(matrix);
            m_State->analyzedSize = matrix.rows();
            m_State->analyzedNonZeros = matrix.nonZeros();
        }

        m_State->solver.factorize(matrix);
        if(m_State->solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Linear solver failed to factorize the system matrix.");
        }

        const Eigen::VectorXd rhs =
          Eigen::VectorXd::Map(rhsVector.data(), static_cast<Eigen::Index>(rhsVector.size()));
        const Eigen::VectorXd result = m_State->solver.solve(rhs);
        if(m_State->solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Linear solver failed to solve the system.");
        }

        return {result.data(), result.data() + result.size()};
    }

    std::vector<double> CLinearSolver::solveEigen(const SquareMatrix & t_MatrixA,
                                                  const std::vector<double> & t_VectorB)
    {
        CLinearSolver oneShot;
        return oneShot.solve(t_MatrixA, t_VectorB);
    }

}   // namespace HygroThermFEM
