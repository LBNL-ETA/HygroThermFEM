#include <iostream>
#include <iomanip>

#include "PrintResults.hxx"

void Helper::printResults(const std::vector<std::vector<double>> & values,
                          std::string_view tableName)
{
    std::cout << "const std::vector<std::vector<double>> " << tableName.data() << " {\n";
    for(auto i = 0u; i < values.size(); ++i)
    {
        std::cout << "   {";
        for(auto j = 0u; j < values[i].size(); ++j)
        {
            std::cout << std::fixed << std::setprecision(6) << values[i][j];
            if(j != values[i].size() - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << "}";
        if(i != values.size() - 1)
        {
            std::cout << ",";
        }
        std::cout << '\n';
    }
    std::cout << "};\n";
}

void Helper::printResults(const std::vector<double> & values, std::string_view tableName)
{
    std::cout << "const std::vector<double> " << tableName.data() << " {";
    for(auto i = 0u; i < values.size(); ++i)
    {
        std::cout << std::scientific << std::setprecision(6) << values[i];
        if(i != values.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "};\n";
}
