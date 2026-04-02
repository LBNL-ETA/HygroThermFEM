#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace TestHelper
{
    inline void printVector(const std::string & name,
                            const std::vector<double> & vec,
                            std::ostream & out = std::cerr)
    {
        out << "    const std::vector<double> " << name << "{";
        for(size_t idx = 0; idx < vec.size(); ++idx)
        {
            if(idx > 0)
            {
                out << ", ";
            }
            // Use scientific for very small or very large values
            if(vec[idx] != 0.0 && (std::abs(vec[idx]) < 1e-4 || std::abs(vec[idx]) >= 1e6))
            {
                out << std::scientific << std::setprecision(6) << vec[idx];
            }
            else
            {
                out << std::fixed << std::setprecision(6) << vec[idx];
            }
        }
        out << "};\n";
    }

    inline void printVector2D(const std::string & name,
                              const std::vector<std::vector<double>> & vec,
                              std::ostream & out = std::cerr)
    {
        out << std::fixed << std::setprecision(6);
        out << "    const std::vector<std::vector<double>> " << name << "{\n";
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
