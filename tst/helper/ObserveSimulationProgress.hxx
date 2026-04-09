#pragma once

#include <vector>

#include "../../src/TimestepObserver.hxx"

namespace TestHelper
{
    //! Counts how many times each non-zero subdivision level fires.
    //! Index 0 corresponds to division level 1, index 1 to level 2, etc.
    class ObserveSimulationProgress : public Timesteps::TimestepObserver
    {
    public:
        void levelChanged(unsigned divisionLevel, unsigned) override
        {
            if(divisionLevel == 0)
            {
                return;
            }
            const size_t idx = divisionLevel - 1;
            if(idx >= m_SimulationCalls.size())
            {
                m_SimulationCalls.resize(idx + 1, 0u);
            }
            ++m_SimulationCalls[idx];
        }

        [[nodiscard]] const std::vector<unsigned> & calls() const
        {
            return m_SimulationCalls;
        }

    private:
        std::vector<unsigned> m_SimulationCalls;
    };
}   // namespace TestHelper
