#pragma once

#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace TestHelper
{
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
