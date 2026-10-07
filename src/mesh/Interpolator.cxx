#include <cmath>
#include <functional>

#include "Interpolator.hxx"

namespace FenestrationCommon
{
    double Interpolator::f(const point &t_point1,
                           const point &t_point2,
                           double t_position)
    {
        const auto delta = t_point2.x - t_point1.x;
        return delta != 0 ? (t_position - t_point1.x) / (t_point2.x - t_point1.x) : 0;
    }

    Interpolator::Interpolator(const Interpolation t_interpolation) :
        m_Interpolation(t_interpolation)
    {}

    double Interpolator::interpolate(const point &t_point1,
                                     const point &t_point2,
                                     double t_position) const
    {
        const double f_pos = f(t_point1, t_point2, t_position);

        switch(m_Interpolation)
        {
            case Interpolation::Linear:
                return f_pos * t_point2.y + (1 - f_pos) * t_point1.y;
            case Interpolation::Logarithmic:
                return std::pow(t_point2.y, f_pos) * std::pow(t_point1.y, (1 - f_pos));
        }

        return 0;
    }
}   // namespace FenestrationCommon
