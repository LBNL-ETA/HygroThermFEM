#include <cmath>
#include <functional>

#include "Interpolator.hxx"

namespace FenestrationCommon
{
    double Interpolator::f(const std::pair<double, double> & t_point1,
                           const std::pair<double, double> & t_point2,
                           const double t_position)
    {
        const auto delta = t_point2.first - t_point1.first;
        return delta != 0 ? (t_position - t_point1.first) / (t_point2.first - t_point1.first) : 0;
    }

    Interpolator::Interpolator(const Interpolation t_interpolation) :
        m_Interpolation(t_interpolation)
    {
        m_Functions[Interpolation::Linear] = [&](const std::pair<double, double> & t_point1,
                                                 const std::pair<double, double> & t_point2,
                                                 const double t_position) {
            const double f_pos = f(t_point1, t_point2, t_position);
            return f_pos * t_point2.second + (1 - f_pos) * t_point1.second;
        };

        m_Functions[Interpolation::Logarithmic] = [&](const std::pair<double, double> & t_point1,
                                                      const std::pair<double, double> & t_point2,
                                                      const double t_position) {
            const double f_pos = f(t_point1, t_point2, t_position);
            return std::pow(t_point2.second, f_pos) * std::pow(t_point1.second, (1 - f_pos));
        };
    }

    double Interpolator::interpolate(const std::pair<double, double> & t_point1,
                                     const std::pair<double, double> & t_point2,
                                     const double t_position) const
    {
        return m_Functions.at(m_Interpolation)(t_point1, t_point2, t_position);
    }
}   // namespace FenestrationCommon
