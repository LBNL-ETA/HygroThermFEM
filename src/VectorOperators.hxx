#pragma once

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

namespace MoisThermFEM
{
    template<class T>
    std::vector<T> operator+(const std::vector<T> & lhs, const std::vector<T> & rhs)
    {
        if(rhs.size() != lhs.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(rhs.size(), 0.0);
        std::transform(rhs.begin(), rhs.end(), lhs.begin(), result.begin(), std::plus<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator-(const std::vector<T> & lhs, const std::vector<T> & rhs)
    {
        if(rhs.size() != lhs.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(rhs.size(), 0);
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), result.begin(), std::minus<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator*(const std::vector<T> & lhs, const std::vector<T> & rhs)
    {
        if(rhs.size() != lhs.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(rhs.size(), 0);
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), result.begin(), std::multiplies<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator*(const std::vector<T> & lhs, const T value)
    {
        std::vector<T> result(lhs.size(), 0);
        std::transform(
          lhs.begin(), lhs.end(), result.begin(), std::bind1st(std::multiplies<T>(), value));

        return result;
    }

    template<class T>
    std::vector<T> operator*(const T value, const std::vector<T> & lhs)
    {
        return lhs * value;
    }

    template<class T>
    std::vector<T> operator/(const std::vector<T> & lhs, const std::vector<T> & rhs)
    {
        if(rhs.size() != lhs.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(rhs.size(), 0);
        std::transform(lhs.begin(), lhs.end(), rhs.begin(), result.begin(), std::divides<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator/(const std::vector<T> & lhs, const T value)
    {
        return lhs * (1 / value);
    }

    template<class T>
    std::vector<T> operator/(const T value, const std::vector<T> & lhs)
    {
        std::vector<T> result(lhs.size(), 0);
        std::transform(lhs.begin(), lhs.end(), result.begin(), [&](T t) { return value / t; });

        return result;
    }

}   // namespace MoisThermFEM