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

        point & operator+(const point & rhs)
        {
            x += rhs.x;
            y += rhs.y;
            return *this;
        }

        point & operator+=(const point & rhs)
        {
            return *this + rhs;
        }

        double x;
        double y;
    };
}   // namespace FenestrationCommon