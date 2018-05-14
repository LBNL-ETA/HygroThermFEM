#ifndef TUPLEASMAPKEY_VARIADICMAP_H
#define TUPLEASMAPKEY_VARIADICMAP_H

#include <map>
#include <Eigen/Sparse>

namespace FenestrationCommon
{
    class DKEigenMap
    {
    public:
        DKEigenMap() = default;

        double & operator()(const size_t & v1, const size_t & v2)
        {
            return m_Types[std::make_pair(v1, v2)];
        }

        double at(const size_t & v1, const size_t & v2) const
        {
            double value{0};
            if(exists(v1, v2))
            {
                value = m_Types.at(std::make_pair(v1, v2));
            }
            return value;
        }

        bool exists(const size_t & v1, const size_t & v2) const
        {
            return m_Types.count(std::make_pair(v1, v2)) != 0;
        }

        void clear()
        {
            m_Types.clear();
        }

        size_t size() const {
            return m_Types.size();
        }

        std::vector<Eigen::Triplet<double>> triplets() const
        {
            std::vector<Eigen::Triplet<double>> triplets;

            for(auto & val : m_Types)
            {
                triplets.emplace_back(int(val.first.first), int(val.first.second), val.second);
            }

            return triplets;
        }

    private:
        std::map<std::pair<size_t, size_t>, double> m_Types;
    };

}   // namespace FenestrationCommon


#endif   // TUPLEASMAPKEY_VARIADICMAP_H
