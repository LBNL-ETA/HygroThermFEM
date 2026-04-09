#pragma once

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace TestHelper
{
    //! Collects violations of physical invariants across a transient simulation.
    //! Each check appends to an internal log; call `report()` to get a summary
    //! and `passed()` to check if all invariants held.
    class PhysicalInvariants
    {
    public:
        //! Every value at every timestep must be >= minVal.
        void expectBoundedBelow(const std::vector<std::vector<double>> & data,
                                const double minVal,
                                const std::string & label)
        {
            for(size_t step = 0; step < data.size(); ++step)
            {
                for(size_t node = 0; node < data[step].size(); ++node)
                {
                    if(data[step][node] < minVal - 1e-10)
                    {
                        addViolation(label, step, node,
                                     "below minimum " + std::to_string(minVal),
                                     data[step][node]);
                    }
                }
            }
        }

        //! Every value at every timestep must be <= maxVal.
        void expectBoundedAbove(const std::vector<std::vector<double>> & data,
                                const double maxVal,
                                const std::string & label)
        {
            for(size_t step = 0; step < data.size(); ++step)
            {
                for(size_t node = 0; node < data[step].size(); ++node)
                {
                    if(data[step][node] > maxVal + 1e-10)
                    {
                        addViolation(label, step, node,
                                     "above maximum " + std::to_string(maxVal),
                                     data[step][node]);
                    }
                }
            }
        }

        //! Every value must stay within [minVal, maxVal].
        void expectBounded(const std::vector<std::vector<double>> & data,
                           const double minVal,
                           const double maxVal,
                           const std::string & label)
        {
            expectBoundedBelow(data, minVal, label);
            expectBoundedAbove(data, maxVal, label);
        }

        //! At each node, value should not decrease from one timestep to the next.
        //! A small tolerance accounts for solver noise.
        void expectNonDecreasing(const std::vector<std::vector<double>> & data,
                                 const double tolerance,
                                 const std::string & label)
        {
            for(size_t step = 1; step < data.size(); ++step)
            {
                for(size_t node = 0; node < data[step].size(); ++node)
                {
                    if(data[step][node] < data[step - 1][node] - tolerance)
                    {
                        addViolation(label, step, node,
                                     "decreased from " + std::to_string(data[step - 1][node]),
                                     data[step][node]);
                    }
                }
            }
        }

        //! At each node, value should not increase from one timestep to the next.
        void expectNonIncreasing(const std::vector<std::vector<double>> & data,
                                 const double tolerance,
                                 const std::string & label)
        {
            for(size_t step = 1; step < data.size(); ++step)
            {
                for(size_t node = 0; node < data[step].size(); ++node)
                {
                    if(data[step][node] > data[step - 1][node] + tolerance)
                    {
                        addViolation(label, step, node,
                                     "increased from " + std::to_string(data[step - 1][node]),
                                     data[step][node]);
                    }
                }
            }
        }

        //! No value should jump by more than maxFraction of its magnitude
        //! between consecutive timesteps.
        void expectSmooth(const std::vector<std::vector<double>> & data,
                          const double maxFraction,
                          const std::string & label)
        {
            for(size_t step = 1; step < data.size(); ++step)
            {
                for(size_t node = 0; node < data[step].size(); ++node)
                {
                    const double prev = data[step - 1][node];
                    const double curr = data[step][node];
                    const double scale = std::max(std::abs(prev), std::abs(curr));
                    if(scale > 1e-10 && std::abs(curr - prev) > maxFraction * scale)
                    {
                        addViolation(label, step, node,
                                     "jump of " + std::to_string(std::abs(curr - prev))
                                       + " (max allowed " + std::to_string(maxFraction * scale) + ")",
                                     curr);
                    }
                }
            }
        }

        //! Symmetric node pairs should have equal values at every timestep.
        //! Pairs are given as indices into the per-step vectors.
        void expectSymmetric(const std::vector<std::vector<double>> & data,
                             const std::vector<std::pair<size_t, size_t>> & pairs,
                             const double tolerance,
                             const std::string & label)
        {
            for(size_t step = 0; step < data.size(); ++step)
            {
                for(const auto & [idx1, idx2] : pairs)
                {
                    if(std::abs(data[step][idx1] - data[step][idx2]) > tolerance)
                    {
                        std::ostringstream msg;
                        msg << "asymmetry: node " << idx1 << " = " << data[step][idx1]
                            << " vs node " << idx2 << " = " << data[step][idx2];
                        m_Violations.push_back(
                          label + " [step " + std::to_string(step) + "]: " + msg.str());
                    }
                }
            }
        }

        [[nodiscard]] bool passed() const
        {
            return m_Violations.empty();
        }

        [[nodiscard]] size_t violationCount() const
        {
            return m_Violations.size();
        }

        [[nodiscard]] std::string report() const
        {
            if(m_Violations.empty())
            {
                return "All physical invariants satisfied.";
            }
            std::ostringstream out;
            out << m_Violations.size() << " violation(s):\n";
            const size_t maxShow = 20;
            for(size_t idx = 0; idx < std::min(m_Violations.size(), maxShow); ++idx)
            {
                out << "  - " << m_Violations[idx] << "\n";
            }
            if(m_Violations.size() > maxShow)
            {
                out << "  ... and " << (m_Violations.size() - maxShow) << " more\n";
            }
            return out.str();
        }

    private:
        void addViolation(const std::string & label,
                          const size_t step,
                          const size_t node,
                          const std::string & detail,
                          const double value)
        {
            m_Violations.push_back(
              label + " [step " + std::to_string(step) + ", node " + std::to_string(node)
              + "]: " + detail + " (value = " + std::to_string(value) + ")");
        }

        std::vector<std::string> m_Violations;
    };
}   // namespace TestHelper
