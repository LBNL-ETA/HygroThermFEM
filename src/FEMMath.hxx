#pragma once

#include <vector>
#include <cmath>

namespace HygroThermFEM
{
    static const double HTFEM_PI = 4.0 * std::atan(1.0);

    template<typename T>
    T norm(const std::vector<T> & t_vector)
    {
        double result{0};
        std::for_each(t_vector.begin(), t_vector.end(), [&](T n) { result += n * n; });

        return std::pow(result, 0.5);
    }

    template<typename T>
    T errorNorm(const std::vector<T> & vec1, const std::vector<T> & vec2)
    {
        auto norm1 = norm(vec1);
        auto norm2 = norm(vec2);
        if(norm1 == 0)
        {
            norm1 = 1e-10;
            if(norm2 == 0)
            {
                norm2 = norm1;
            }
        }

        return std::abs(norm1 - norm2) / norm1;
    }

    double radians(double d);

    double degrees(double r);
}   // namespace HygroThermFEM
