#ifndef VECTOR_OPERATORS_HXX
#define VECTOR_OPERATORS_HXX

#include <algorithm>

namespace MoisThermFEM
{
    /// Operator +

    std::vector<double> operator+(const std::vector<double> & first,
                                  const std::vector<double> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<double> result(first.size(), 0);
        std::transform(first.begin(), first.end(), second.begin(), result.begin(), std::plus<double>());

        return result;
    }

    std::vector<double> operator+(const std::vector<double> & first, const double second)
    {
        std::vector<double> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), std::bind1st(std::plus<double>(), second));

        return result;
    }

    std::vector<double> operator+(const double first, const std::vector<double> & second)
    {
        return operator+(second, first);
    }

    /// Operator -

    std::vector<double> operator-(const std::vector<double> & first,
                                  const std::vector<double> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<double> result(first.size(), 0);
        std::transform(first.begin(), first.end(), second.begin(), result.begin(), std::minus<double>());

        return result;
    }

    std::vector<double> operator-(const std::vector<double> & first, const double second)
    {
        std::vector<double> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), std::bind1st(std::minus<double>(), second));

        return result;
    }

    std::vector<double> operator-(const double first, const std::vector<double> & second)
    {
        return operator-(second, first);
    }

    /// Operator *

    std::vector<double> operator*(const std::vector<double> & first,
                                  const std::vector<double> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<double> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), second.begin(), result.begin(), std::multiplies<double>());

        return result;
    }

    std::vector<double> operator*(const std::vector<double> & first, const double second)
    {
        std::vector<double> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), std::bind1st(std::multiplies<double>(), second));

        return result;
    }

    std::vector<double> operator*(const double first, const std::vector<double> & second)
    {
        return operator*(second, first);
    }

    /// Operator /

    std::vector<double> operator/(const std::vector<double> & first,
                                  const std::vector<double> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<double> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), second.begin(), result.begin(), std::divides<double>());

        return result;
    }

    std::vector<double> operator/(const std::vector<double> & first, const double second)
    {
        std::vector<double> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), std::bind1st(std::divides<double>(), second));

        return result;
    }

    std::vector<double> operator/(const double first, const std::vector<double> & second)
    {
        std::vector<double> result(second.size(), 0);
        std::transform(
          second.begin(), second.end(), result.begin(), [&](double t) { return first / t; });

        return result;
    }
}   // namespace MoisThermFEM

#endif