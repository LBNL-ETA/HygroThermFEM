#include "FEMMath.hxx"

#include <cmath>
#include <algorithm>

namespace HygroThermFEM
{
    double norm(const std::vector<double> & t_vector)
    {
        double result{0};
        std::for_each(t_vector.begin(), t_vector.end(), [&](double n) { result += n * n; });

        return std::pow(result, 0.5);
    }

    double radians(const double d)
    {
        return d * HTFEM_PI / 180;
    }

    double degrees(const double r)
    {
        return r * 180 / HTFEM_PI;
    }
}   // namespace HygroThermFEM
