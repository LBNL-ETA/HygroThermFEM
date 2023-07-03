#include "FEMMath.hxx"

#include <cmath>
#include <algorithm>

namespace HygroThermFEM
{
    double radians(const double d)
    {
        return d * HTFEM_PI / 180;
    }

    double degrees(const double r)
    {
        return r * 180 / HTFEM_PI;
    }
}   // namespace HygroThermFEM
