#pragma once

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

namespace HygroThermFEM
{
    /// Operator +

    template<class T>
    std::vector<T> operator+(const std::vector<T> & first, const std::vector<T> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(first.size(), 0);
        std::transform(first.begin(), first.end(), second.begin(), result.begin(), std::plus<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator+(const std::vector<T> & first, const T second)
    {
        std::vector<T> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), std::bind(std::plus<T>(), second));

        return result;
    }

    template<class T>
    std::vector<T> operator+(const T first, const std::vector<T> & second)
    {
        return operator+(second, first);
    }

    /// Operator -
    template<class T>
    std::vector<T> operator-(const std::vector<T> & first, const std::vector<T> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(first.size(), 0);
        std::transform(first.begin(), first.end(), second.begin(), result.begin(), std::minus<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator-(const std::vector<T> & first, const T second)
    {
        std::vector<T> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), std::bind(std::minus<T>(), second));

        return result;
    }

    template<class T>
    std::vector<T> operator-(const T first, const std::vector<T> & second)
    {
        return operator-(second, first);
    }

    /// Operator *

    template<class T>
    std::vector<T> operator*(const std::vector<T> & first, const std::vector<T> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), second.begin(), result.begin(), std::multiplies<T>());

        return result;
    }

    template<class T>
    std::vector<T> operator*(const std::vector<T> & first, const T second)
    {
        std::vector<T> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), result.begin(), [&](auto && data) { return data * second; });

        return result;
    }

    template<class T>
    std::vector<T> operator*(const T first, const std::vector<T> & second)
    {
        return operator*(second, first);
    }

    /// Operator /

    template<class T>
    std::vector<T> operator/(const std::vector<T> & first, const std::vector<T> & second)
    {
        if(first.size() != second.size())
        {
            throw std::runtime_error("Vectors must be identical in size.");
        }

        std::vector<T> result(first.size(), 0);
        std::transform(
          first.begin(), first.end(), second.begin(), result.begin(), std::divides<T>());

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

}   // namespace HygroThermFEM
