#pragma once

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace TestHelper
{
    //! Opt-in CSV capture for the Python reference-solver comparison (hygrothermfem_python).
    //!
    //! A complete no-op unless the HTF_DUMP_DIR environment variable is set, so CI
    //! and normal test runs write nothing. Rows accumulate in memory and the file is
    //! written on destruction, which keeps every call site to one line and never
    //! leaves a partial file behind a failed test.
    //!
    //! Emits exactly the two layouts the Python loader reads (refdata.py _READERS):
    //!   step_node        header "step,x0,x1,..."     one row per timestep
    //!   step_field_node  header "step,field,n0,..."  one row per (timestep, field)
    //! The layout is fixed by the first addRow call; mixing the two in one dump
    //! throws, as does a row whose value count differs from the declared column
    //! count. Values are written with 12 significant digits: the comparisons
    //! resolve to 1e-5 and the stored tolerances are 5e-5, so anything fewer would
    //! put rounding into the same range as the result.
    class CsvDump
    {
    public:
        CsvDump(const std::string & fileName, std::size_t nColumns) : m_NColumns(nColumns)
        {
            if(const char * dir = std::getenv("HTF_DUMP_DIR"); dir != nullptr && *dir != '\0')
            {
                m_Path = std::filesystem::path(dir) / fileName;
            }
        }

        CsvDump(const CsvDump &) = delete;
        CsvDump & operator=(const CsvDump &) = delete;

        //! step_node layout: one row per timestep.
        void addRow(std::size_t step, const std::vector<double> & values)
        {
            appendRow(step, "", values, false);
        }

        //! step_field_node layout: one row per (timestep, field).
        void addRow(std::size_t step, const std::string & field, const std::vector<double> & values)
        {
            appendRow(step, field, values, true);
        }

        ~CsvDump()
        {
            if(m_Path.empty() || m_Rows.empty())
            {
                return;
            }
            try
            {
                writeFile();
            }
            catch(...)
            {
                // A dump helper must never turn a passing test into a crash; a
                // missing file is diagnosis enough that the write failed.
            }
        }

    private:
        struct Row
        {
            std::size_t step;
            std::string field;
            std::vector<double> values;
        };

        void appendRow(std::size_t step,
                       const std::string & field,
                       const std::vector<double> & values,
                       bool fieldLayout)
        {
            if(m_Path.empty())
            {
                return;
            }
            if(values.size() != m_NColumns)
            {
                throw std::invalid_argument("CsvDump: row has " + std::to_string(values.size())
                                            + " values, expected " + std::to_string(m_NColumns));
            }
            if(!m_Rows.empty() && fieldLayout != m_FieldLayout)
            {
                throw std::logic_error("CsvDump: cannot mix step_node and step_field_node rows");
            }
            m_FieldLayout = fieldLayout;
            m_Rows.push_back({step, field, values});
        }

        void writeFile() const
        {
            std::filesystem::create_directories(m_Path.parent_path());
            std::ofstream out(m_Path);
            out << std::setprecision(12);
            out << (m_FieldLayout ? "step,field" : "step");
            for(std::size_t col = 0; col < m_NColumns; ++col)
            {
                out << ",x" << col;
            }
            out << "\n";
            for(const auto & row : m_Rows)
            {
                out << row.step;
                if(m_FieldLayout)
                {
                    out << "," << row.field;
                }
                for(const auto & value : row.values)
                {
                    out << "," << value;
                }
                out << "\n";
            }
        }

        std::filesystem::path m_Path;   // empty => dumping disabled
        std::size_t m_NColumns{0};
        bool m_FieldLayout{false};
        std::vector<Row> m_Rows;
    };

    //! Bottom-node-row extraction for the one-element-high beam/slab 1D proxy.
    //!
    //! Nodes are 1-based and column-major with nodesPerColumn nodes per x-column
    //! (SlabBuilder/BeamBuilder order), so the bottom node of column c sits at
    //! index c * nodesPerColumn in the 0-based solution vector.
    inline std::vector<double> bottomRow(const std::vector<double> & nodal,
                                         std::size_t nColumns,
                                         std::size_t nodesPerColumn)
    {
        std::vector<double> row;
        row.reserve(nColumns);
        for(std::size_t col = 0; col < nColumns; ++col)
        {
            row.push_back(nodal[col * nodesPerColumn]);
        }
        return row;
    }
}   // namespace TestHelper
