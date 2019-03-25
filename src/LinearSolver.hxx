#pragma once

#include <vector>

#include "SquareMatrix.hxx"

namespace HygroThermFEM
{
    class CLinearSolver
    {
    public:
        CLinearSolver() = default;

        static std::vector<double> solveEigen(const SquareMatrix & t_MatrixA,
                                              const std::vector<double> & t_VectorB);

    private:
        std::vector<double> checkSingularity(SquareMatrix & t_MatrixA) const;
        std::vector<std::vector<size_t>> m_RowIndexes;
        std::vector<std::vector<size_t>> m_ColumnIndexes;
    };
}   // namespace HygroThermFEM
