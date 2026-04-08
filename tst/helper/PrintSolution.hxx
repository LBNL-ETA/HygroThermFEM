#pragma once
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace TestHelper
{
    namespace detail
    {
        //! Returns the C++ type name string used when emitting a vector
        //! literal for type T. Specialise for any new type that should be
        //! supported by printVector / printVector2D.
        template<typename T>
        constexpr const char * typeName() = delete;

        template<>
        constexpr const char * typeName<double>()
        {
            return "double";
        }
        template<>
        constexpr const char * typeName<int>()
        {
            return "int";
        }
        template<>
        constexpr const char * typeName<unsigned>()
        {
            return "unsigned";
        }
        template<>
        constexpr const char * typeName<std::size_t>()
        {
            return "size_t";
        }

        //! Writes a single scalar value to `out`, using scientific format
        //! only for floating-point values that are very small or very large.
        template<typename T>
        void writeScalar(std::ostream & out, const T & value)
        {
            if constexpr(std::is_floating_point_v<T>)
            {
                if(value != T{0} && (std::abs(value) < 1e-4 || std::abs(value) >= 1e6))
                {
                    out << std::scientific << std::setprecision(6) << value;
                }
                else
                {
                    out << std::fixed << std::setprecision(6) << value;
                }
            }
            else
            {
                out << value;
            }
        }
    }   // namespace detail

    template<typename T>
    void printVector(const std::string & name,
                     const std::vector<T> & vec,
                     std::ostream & out = std::cerr)
    {
        out << "    const std::vector<" << detail::typeName<T>() << "> " << name << "{";
        for(size_t idx = 0; idx < vec.size(); ++idx)
        {
            if(idx > 0)
            {
                out << ", ";
            }
            detail::writeScalar(out, vec[idx]);
        }
        out << "};\n";
    }

    template<typename T>
    void printVector2D(const std::string & name,
                       const std::vector<std::vector<T>> & vec,
                       std::ostream & out = std::cerr)
    {
        if constexpr(std::is_floating_point_v<T>)
        {
            out << std::fixed << std::setprecision(6);
        }
        out << "    const std::vector<std::vector<" << detail::typeName<T>() << ">> " << name
            << "{\n";
        for(size_t idx = 0; idx < vec.size(); ++idx)
        {
            out << "      {";
            for(size_t jdx = 0; jdx < vec[idx].size(); ++jdx)
            {
                if(jdx > 0)
                {
                    out << ", ";
                }
                out << vec[idx][jdx];
            }
            out << "}";
            if(idx + 1 < vec.size())
            {
                out << ",";
            }
            out << "\n";
        }
        out << "    };\n";
    }

    inline void printUnsigned(const std::string & name,
                              const unsigned val,
                              std::ostream & out = std::cerr)
    {
        out << "    " << name << " = " << val << "u;\n";
    }

    //! Opens a file and dumps all solution data in one call.
    //! Usage:
    //!   TestHelper::dumpToFile("results.txt",
    //!       {"waterContent", waterContentSolution},
    //!       {"temperature", temperatureSolution});
    struct NamedVector2D
    {
        std::string name;
        const std::vector<std::vector<double>> & data;
    };

    inline void dumpToFile(const std::string & filename,
                           std::initializer_list<NamedVector2D> vectors)
    {
        std::ofstream out(filename);
        for(const auto & [name, data] : vectors)
        {
            printVector2D(name, data, out);
            out << "\n";
        }
    }
}   // namespace TestHelper
