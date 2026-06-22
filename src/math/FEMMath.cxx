#include "FEMMath.hxx"

#include <cmath>
#include <algorithm>
#include <numbers>

namespace HygroThermFEM
{
    double norm(const std::vector<double> & t_vector)
    {
        double result{0};
        std::ranges::for_each(t_vector, [&](const double n) { result += n * n; });

        return std::pow(result, 0.5);
    }

    double radians(const double d)
    {
        return d * std::numbers::pi / 180;
    }
}   // namespace HygroThermFEM
