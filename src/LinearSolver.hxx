#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>

#include "SquareMatrix.hxx"

namespace FenestrationCommon
{
    class CLinearSolver
    {
    public:
        CLinearSolver() = default;

        static std::vector<double> solveEigen( const SquareMatrix & t_MatrixA,
											   const std::vector< double > & t_VectorB );

    private:
        std::vector<double> checkSingularity(SquareMatrix & t_MatrixA) const;
        std::vector<std::vector<size_t>> m_RowIndexes;
        std::vector<std::vector<size_t>> m_ColumnIndexes;
    };
}   // namespace FenestrationCommon

#endif
