#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include <vector>

#include "SparceSquareMatrix.hxx"

namespace FenestrationCommon
{
    class CLinearSolver
    {
    public:
        CLinearSolver() = default;

        static std::vector<double> solveSystem(SparceSquareMatrix<double> t_MatrixA,
                                               std::vector<double> & t_VectorB);

        static std::vector<double> solveEigenSparse( const SparceSquareMatrix< double > & t_MatrixA,
													 const std::vector< double > & t_VectorB );

    private:
        std::vector<double> checkSingularity(SparceSquareMatrix<double> & t_MatrixA) const;
        std::vector<std::vector<size_t>> m_RowIndexes;
        std::vector<std::vector<size_t>> m_ColumnIndexes;
    };
}   // namespace FenestrationCommon

#endif
