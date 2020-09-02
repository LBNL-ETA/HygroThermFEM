#pragma once

#include <vector>
#include <cmath>

namespace HygroThermFEM
{
    static const double HTFEM_PI = 4.0 * std::atan(1.0);

    double norm(const std::vector<double> & t_vector);

    double radians(double d);

    double degrees(double r);
}   // namespace HygroThermFEM
