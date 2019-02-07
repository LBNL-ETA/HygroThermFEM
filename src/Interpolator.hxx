#pragma once

#include <functional>
#include <map>
#include <utility>

namespace FenestrationCommon {

enum class Interpolation { Linear, Logarithmic };

class Interpolator {
public:
  Interpolator(Interpolation t_interpolation);

  double interpolate(const std::pair<double, double> &t_point1,
                     const std::pair<double, double> &t_point2,
                     double t_position) const;

private:
  static double f(const std::pair<double, double> &t_point1,
           const std::pair<double, double> &t_point2,
           double t_position);

  std::map<Interpolation,
           std::function<double(const std::pair<double, double> &t_point1,
                                const std::pair<double, double> &t_point2,
                                double t_position)>>
      m_Functions;

  Interpolation m_Interpolation;
};

} // namespace FenestrationCommon