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

    void Settings::setMaxDivisions(const unsigned mMaxDivisions)
    {
        m_MaxDivisions = mMaxDivisions;
    }

    unsigned Settings::getNumberOfSubtimesteps() const
    {
        return m_NumberOfSubsegments;
    }

    void Settings::setNumberOfSubtimesteps(const unsigned numberOfSubtimesteps)
    {
        m_NumberOfSubsegments = numberOfSubtimesteps;
    }
}   // namespace Timesteps
