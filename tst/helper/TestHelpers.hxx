#pragma once

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "HygroThermFEM2D.hxx"

namespace TestHelper
{
    //! Golden-value regeneration support. When the environment variable HTF_DUMP_GOLDEN is
    //! set, each dumpGolden call prints the ACTUAL data as a C++ initializer list tagged
    //! with the golden variable's name, so a re-baselining script can patch the hardcoded
    //! expected arrays in the test sources. With the variable unset the calls are no-ops,
    //! so they can stay in the tests permanently.
    inline bool goldenDumpEnabled()
    {
        return std::getenv("HTF_DUMP_GOLDEN") != nullptr;
    }

    inline void dumpGolden(const std::string & name, const std::vector<double> & values)
    {
        if(!goldenDumpEnabled())
        {
            return;
        }
        std::cout << "[GOLDEN] " << name << "\n" << std::setprecision(12) << "{";
        for(size_t idx = 0; idx < values.size(); ++idx)
        {
            std::cout << (idx == 0 ? "" : ", ") << values[idx];
        }
        std::cout << "}\n[/GOLDEN]\n";
    }

    inline void dumpGolden(const std::string & name,
                           const std::vector<std::vector<double>> & values)
    {
        if(!goldenDumpEnabled())
        {
            return;
        }
        std::cout << "[GOLDEN] " << name << "\n" << std::setprecision(12) << "{";
        for(size_t row = 0; row < values.size(); ++row)
        {
            std::cout << (row == 0 ? "{" : ",\n {");
            for(size_t col = 0; col < values[row].size(); ++col)
            {
                std::cout << (col == 0 ? "" : ", ") << values[row][col];
            }
            std::cout << "}";
        }
        std::cout << "}\n[/GOLDEN]\n";
    }

    inline void dumpGolden(const std::string & name,
                           const std::vector<std::vector<HygroThermFEM::NodeFlux>> & values)
    {
        if(!goldenDumpEnabled())
        {
            return;
        }
        std::cout << "[GOLDEN] " << name << "\n" << std::setprecision(12) << "{";
        for(size_t row = 0; row < values.size(); ++row)
        {
            std::cout << (row == 0 ? "{" : ",\n {");
            for(size_t col = 0; col < values[row].size(); ++col)
            {
                std::cout << (col == 0 ? "" : ", ") << "{" << values[row][col].x << ", "
                          << values[row][col].y << "}";
            }
            std::cout << "}";
        }
        std::cout << "}\n[/GOLDEN]\n";
    }

    inline void expectNear(const std::vector<std::vector<double>> & expected,
                             const std::vector<std::vector<double>> & actual,
                             const double tolerance)
    {
        ASSERT_EQ(expected.size(), actual.size());
        for(size_t row = 0; row < expected.size(); ++row)
        {
            ASSERT_EQ(expected[row].size(), actual[row].size()) << "Row " << row;
            for(size_t col = 0; col < expected[row].size(); ++col)
            {
                EXPECT_NEAR(expected[row][col], actual[row][col], tolerance)
                  << "at [" << row << "][" << col << "]";
            }
        }
    }

    inline void expectNear(const std::vector<double> & expected,
                             const std::vector<double> & actual,
                             const double tolerance)
    {
        ASSERT_EQ(expected.size(), actual.size());
        for(size_t idx = 0; idx < expected.size(); ++idx)
        {
            EXPECT_NEAR(expected[idx], actual[idx], tolerance) << "at [" << idx << "]";
        }
    }
}   // namespace TestHelper
