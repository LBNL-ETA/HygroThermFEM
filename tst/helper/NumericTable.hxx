#pragma once

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace TestHelper
{
    //! Rows of numbers read from a delimited text file. Commas and whitespace both
    //! separate fields, `headerLines` leading lines are skipped, blank lines are ignored,
    //! and a field that does not parse as a number reads as NaN so that a sparse
    //! spreadsheet export keeps its column positions.
    using NumericTable = std::vector<std::vector<double>>;

    [[nodiscard]] inline std::vector<double> parseNumericLine(const std::string & line)
    {
        std::string spaced{line};
        std::ranges::replace(spaced, ',', ' ');
        std::istringstream stream{spaced};
        std::vector<double> row;
        std::string field;
        while(stream >> field)
        {
            try
            {
                row.push_back(std::stod(field));
            }
            catch(const std::invalid_argument &)
            {
                row.push_back(std::numeric_limits<double>::quiet_NaN());
            }
        }
        return row;
    }

    [[nodiscard]] inline NumericTable loadNumericTable(const std::string & path,
                                                       const std::size_t headerLines)
    {
        std::ifstream file{path};
        if(!file)
        {
            throw std::runtime_error("Cannot open data file: " + path);
        }
        NumericTable table;
        std::string line;
        std::size_t lineIndex{0u};
        while(std::getline(file, line))
        {
            ++lineIndex;
            if(lineIndex <= headerLines || line.find_first_not_of(" \t\r\n,") == std::string::npos)
            {
                continue;
            }
            table.push_back(parseNumericLine(line));
        }
        return table;
    }

    //! One column of a table.
    [[nodiscard]] inline std::vector<double> column(const NumericTable & table,
                                                    const std::size_t index)
    {
        std::vector<double> values;
        values.reserve(table.size());
        for(const auto & row : table)
        {
            values.push_back(index < row.size() ? row[index]
                                                : std::numeric_limits<double>::quiet_NaN());
        }
        return values;
    }

    //! Linear interpolation of `values` sampled at ascending `coords`, clamped at the ends.
    [[nodiscard]] inline double interpolateAt(const std::vector<double> & coords,
                                              const std::vector<double> & values,
                                              const double position)
    {
        if(position <= coords.front())
        {
            return values.front();
        }
        if(position >= coords.back())
        {
            return values.back();
        }
        const auto upper = std::ranges::upper_bound(coords, position) - coords.begin();
        const auto right = static_cast<std::size_t>(upper);
        const auto left = right - 1u;
        const double fraction{(position - coords[left]) / (coords[right] - coords[left])};
        return values[left] + fraction * (values[right] - values[left]);
    }
}   // namespace TestHelper
