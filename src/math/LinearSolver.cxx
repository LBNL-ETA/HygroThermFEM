#pragma warning(push, 0)
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#pragma warning(pop)

#include <cmath>
#include <fstream>
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

        // TEMP diagnostic for IGU/HygroThermFEM steady crash: is the system non-finite (NaN/Inf from
        // e.g. enclosure radiosity) or singular (a zero-diagonal / unconstrained node)? Remove with
        // the rest of the steady debug instrumentation.
        {
            int nonFiniteA = 0;
            int zeroDiagonal = 0;
            int nonFiniteB = 0;
            for(int col = 0; col < matrix.outerSize(); ++col)
            {
                for(Matrix::InnerIterator it(matrix, col); it; ++it)
                {
                    if(!std::isfinite(it.value()))
                    {
                        ++nonFiniteA;
                    }
                }
            }
            for(int row = 0; row < matrix.rows(); ++row)
            {
                if(matrix.coeff(row, row) == 0.0)
                {
                    ++zeroDiagonal;
                }
            }
            for(const auto value : t_VectorB)
            {
                if(!std::isfinite(value))
                {
                    ++nonFiniteB;
                }
            }
            std::ofstream log("D:\\tmp\\htf_steady.log", std::ios::app);
            log << "    [solve] N=" << matrix.rows() << " nonFiniteA=" << nonFiniteA
                << " nonFiniteB=" << nonFiniteB << " zeroDiagonalRows=" << zeroDiagonal << "\n";
            log.flush();
        }

        Eigen::SparseLU<Matrix> solver;
        solver.analyzePattern(matrix);
        solver.factorize(matrix);

        {
            std::ofstream log("D:\\tmp\\htf_steady.log", std::ios::app);
            log << "    [solve] factorize info=" << static_cast<int>(solver.info()) << "\n";
            log.flush();
        }

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
