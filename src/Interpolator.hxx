#pragma once

#include <functional>
#include <map>
#include <utility>

#include "Point.hxx"

namespace FenestrationCommon
{
    enum class Interpolation
    {
        Linear,
        Logarithmic
    };

    class Interpolator
    {
    public:
        Interpolator(Interpolation t_interpolation);

        double interpolate(const point & t_point1, const point & t_point2, double t_position) const;

    private:
        static double f(const point & t_point1, const point & t_point2, double t_position);

        std::map<
          Interpolation,
          std::function<double(const point & t_point1, const point & t_point2, double t_position)>>
          m_Functions;

        Interpolation m_Interpolation;
    };

}   // namespace FenestrationCommon