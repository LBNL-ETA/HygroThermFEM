#pragma once

namespace Timesteps
{
    //! \brief Singleton class to keep basic settings for timestep division algorithm
    class Settings
    {
    public:
        //! \brief Calling instance of singleton class
        static Settings & Instance();

        //! \brief Retuns maximum number of timestep divisions.
        //!
        //! If engine fails to converge, it will first try to divide current timestep into number of
        //! smaller segments. It will keep doing that till maximum number of divisions has been
        //! reached.
        //! \return Current value of maximum allowed timestep subdivisions.
        [[nodiscard]] unsigned getMaxDivisions() const;

        //! \brief Number of smaller timesteps that program will divide current timestep in case of
        //! convergence failure
        //!
        //! \return Number of timesteps.
        [[nodiscard]] unsigned getNumberOfSubtimesteps() const;

    private:
        Settings() = default;
        ~Settings() = default;

        unsigned m_MaxDivisions{3};
        unsigned m_NumberOfSubsegments{10};
    };
}   // namespace Timesteps