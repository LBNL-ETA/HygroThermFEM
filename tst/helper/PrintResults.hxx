#pragma once

#include <vector>

namespace Helper
{
    void printResults(const std::vector<std::vector<double>> & values, std::string_view tableName);
    void printResults(const std::vector<double> & values, std::string_view tableName);
}   // namespace Helper