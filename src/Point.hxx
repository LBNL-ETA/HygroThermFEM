#pragma once

namespace FenestrationCommon
{
    //! \brief Simple structure used to store points.
    struct point
    {
        point(const double x, const double y) : x(x), y(y){};

        bool operator==(const point & rhs) const
        {
            return x == rhs.x && y == rhs.y;
        }

        bool operator!=(const point & rhs) const
        {
            return !(rhs == *this);
        }

        double x;
        double y;
    };
}   // namespace FenestrationCommon