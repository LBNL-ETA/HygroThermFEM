#include <stdexcept>
#include <cassert>
#include <cmath>

#pragma warning(push, 0)
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#pragma warning(pop) 

#include "LinearSolver.hxx"

using FenestrationCommon::SquareMatrix;

namespace FenestrationCommon
{
    std::vector<double>
      CLinearSolver::checkSingularity(SquareMatrix & t_MatrixA) const
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

    std::vector<double> CLinearSolver::solveEigen( const SquareMatrix & t_MatrixA,
												   const std::vector< double > & t_VectorB )
    {
        using Matrix = Eigen::MatrixXd;
        using Vector = Eigen::VectorXd;

        const auto size = t_MatrixA.size();


        Matrix A(t_MatrixA.m_Matrix);
        Vector B(size);

        for(auto j = 0u; j < t_VectorB.size(); ++j)
        {
            B[j] = t_VectorB[j];
        }

        //Eigen::SparseLU<Matrix> solver;
        //solver.analyzePattern(A);
        //solver.factorize(A);
        //Vector y = solver.solve(B);
		Vector y = A.colPivHouseholderQr().solve(B);

        std::vector<double> solution(y.size());
        for(auto i = 0u; i < size; ++i)
        {
            solution[i] = y[i];
        }

        return solution;
    }

}   // namespace FenestrationCommon
