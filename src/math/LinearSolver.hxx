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
    };
}   // namespace HygroThermFEM
