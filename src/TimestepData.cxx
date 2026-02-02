#include "TimestepData.hxx"

namespace Timesteps
{
    Settings & Settings::Instance()
    {
        static Settings m_Instance;
        return m_Instance;
    }

    unsigned Settings::getMaxDivisions() const
    {
        return m_MaxDivisions;
    }

    unsigned Settings::getNumberOfSubtimesteps() const
    {
        return m_NumberOfSubsegments;
    }
}   // namespace Timesteps
