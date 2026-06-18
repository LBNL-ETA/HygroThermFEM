#include <array>

#include "FeKernels.hxx"

namespace lbnl::errorest::detail
{
    namespace
    {
        //! 1 / sqrt(3) -- the bilinear-quad Gauss point coordinate.
        constexpr double kGaussCoord = 0.5773502691896257;
    }

    std::array<double, 4> shapeFunctions(const double rrr, const double sss)
    {
        return {
            0.25 * (1.0 - rrr) * (1.0 - sss),
            0.25 * (1.0 + rrr) * (1.0 - sss),
            0.25 * (1.0 + rrr) * (1.0 + sss),
            0.25 * (1.0 - rrr) * (1.0 + sss),
        };
    }

    std::array<double, 4> recoveryBasis(const Point & pnt)
    {
        return {1.0, pnt[0], pnt[1], pnt[0] * pnt[1]};
    }

    std::array<Point, 4> gaussRule()
    {
        return {{
            {-kGaussCoord, -kGaussCoord},
            {kGaussCoord, -kGaussCoord},
            {kGaussCoord, kGaussCoord},
            {-kGaussCoord, kGaussCoord},
        }};
    }

    double jacobianDet(const std::array<double, 4> & xxx,
                       const std::array<double, 4> & yyy,
                       const double sss,
                       const double ttt)
    {
        const std::array<std::array<double, 4>, 4> mat = {{
            {0.0, 1.0 - ttt, ttt - sss, sss - 1.0},
            {ttt - 1.0, 0.0, sss + 1.0, -ttt - sss},
            {sss - ttt, -sss - 1.0, 0.0, ttt + 1.0},
            {1.0 - sss, sss + ttt, -1.0 - ttt, 0.0},
        }};

        std::array<double, 4> tmp = {0.0, 0.0, 0.0, 0.0};
        for (std::size_t col = 0; col < 4; ++col)
        {
            for (std::size_t row = 0; row < 4; ++row)
            {
                tmp[col] += xxx[row] * mat[row][col];
            }
        }

        double det = 0.0;
        for (std::size_t idx = 0; idx < 4; ++idx)
        {
            det += tmp[idx] * yyy[idx];
        }
        return det * 0.125;
    }

    Flux applyTensor(const Tensor & ten, const Flux & vec)
    {
        return {
            ten[0] * vec[0] + ten[1] * vec[1],
            ten[2] * vec[0] + ten[3] * vec[1],
        };
    }

    double dot2(const Flux & lhs, const Flux & rhs)
    {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1];
    }

    double dot4(const std::array<double, 4> & lhs, const std::array<double, 4> & rhs)
    {
        double acc = 0.0;
        for (std::size_t idx = 0; idx < 4; ++idx)
        {
            acc += lhs[idx] * rhs[idx];
        }
        return acc;
    }
}
