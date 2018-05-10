#include <stdexcept>
#include <cassert>
#include <cmath>

#pragma warning(push, 0)
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#pragma warning(pop) 

#include "LinearSolver.hxx"

using FenestrationCommon::SparceSquareMatrix;

namespace FenestrationCommon
{
    std::vector<double>
      CLinearSolver::checkSingularity(SparceSquareMatrix<double> & t_MatrixA) const
    {
        const auto size = t_MatrixA.size();
        std::vector<double> vv;

        for(size_t i = 0; i < size; ++i)
        {
            double aamax = 0;
            for(size_t j = 0; j < size; ++j)
            {
                const auto absCellValue = fabs(t_MatrixA(i, j));
                if(absCellValue > aamax)
                {
                    aamax = absCellValue;
                }
            }
            if(aamax == 0)
            {
                assert(aamax != 0);
            }
            vv.push_back(1 / aamax);
        }

        return vv;
    }

    std::vector<double> CLinearSolver::solveSystem(SparceSquareMatrix<double> t_MatrixA,
                                                   std::vector<double> & t_VectorB)
    {
        std::vector<double> solution(t_VectorB.size());

        // Examine zeros on diagonal
        const auto tiny = 1e-9;   // Need to calculate this number
        for(auto i = 0u; i < t_MatrixA.size(); ++i)
        {
            if(t_MatrixA(i, i) == 0)
            {
                t_MatrixA(i, i) = tiny;
            }
        }

        t_MatrixA.makeUpperTriangular(t_VectorB);

        for(auto i = t_MatrixA.size(); i-- > 0;)
        {
            for(auto j = t_MatrixA.size() - 1; j > i; --j)
            {
                solution[i] -= t_MatrixA(i, j) * solution[j];
            }
            solution[i] += t_VectorB[i];
            solution[i] = solution[i] / t_MatrixA(i, i);
        }

        return solution;
    }

    std::vector<double>
      CLinearSolver::solveEigenSparse(const SparceSquareMatrix<double> & t_MatrixA,
                                      const std::vector<double> & t_VectorB)
    {
        using SpMatrix = Eigen::SparseMatrix<double>;
        using Vector = Eigen::VectorXd;

        const auto size = t_MatrixA.size();

        SpMatrix A(size, size);
        Vector B(size);

        auto coefficients = t_MatrixA.triplets();

        A.setFromTriplets(coefficients.begin(), coefficients.end());

        for(auto j = 0u; j < t_VectorB.size(); ++j)
        {
            B[j] = t_VectorB[j];
        }

        Eigen::SparseLU<SpMatrix> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        Vector y = solver.solve(B);

        std::vector<double> solution(y.size());
        for(auto i = 0u; i < size; ++i)
        {
            solution[i] = y[i];
        }

        return solution;
    }

}   // namespace FenestrationCommon
