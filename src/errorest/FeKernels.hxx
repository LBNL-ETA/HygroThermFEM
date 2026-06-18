#pragma once

#include <array>

#include "ErrorEstimator.hxx"

//! Pure finite-element geometry kernels for the bilinear quadrilateral.
//! These are the domain-independent building blocks; they hold no state and are
//! individually unit-testable. The linear algebra (SPR normal-equation solve)
//! lives in the orchestration unit and reuses Windows-CalcEngine matrices.
namespace lbnl::errorest::detail
{
    //! Bilinear quad shape functions N_k at reference point (rrr, sss).
    [[nodiscard]] std::array<double, 4> shapeFunctions(double rrr, double sss);

    //! SPR polynomial basis p = [1, x, y, x*y] evaluated at a global point.
    [[nodiscard]] std::array<double, 4> recoveryBasis(const Point & pnt);

    //! The 2x2 Gauss rule (points +/- 1/sqrt(3), unit weights) on the reference square.
    [[nodiscard]] std::array<Point, 4> gaussRule();

    //! Jacobian determinant of the bilinear iso-parametric map at (sss, ttt).
    [[nodiscard]] double jacobianDet(const std::array<double, 4> & xxx,
                                     const std::array<double, 4> & yyy,
                                     double sss,
                                     double ttt);

    //! Apply a row-major 2x2 tensor to a flux vector.
    [[nodiscard]] Flux applyTensor(const Tensor & ten, const Flux & vec);

    //! Dot product of two flux vectors.
    [[nodiscard]] double dot2(const Flux & lhs, const Flux & rhs);

    //! Dot product of two 4-vectors (basis . coefficients).
    [[nodiscard]] double dot4(const std::array<double, 4> & lhs,
                              const std::array<double, 4> & rhs);
}
