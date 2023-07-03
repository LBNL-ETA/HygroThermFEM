#pragma once

#include "MultiDomain.hxx"
namespace Timesteps
{
    class TimestepObserver;
}

namespace HygroThermFEM
{
    //! \brief Assign observer to thermal part of the engine.
    //!
    //! @param observer observer that will be notified about each timestep
    void subscribeThermal(MultiDomain & domain, Timesteps::TimestepObserver * observer);

    //! \brief Unsubscribe from thermal notifications
    //!
    //! @param observer observer that will be notified about each timestep
    void unsubscribeThermal(MultiDomain & domain, Timesteps::TimestepObserver * observer);

    //! \brief Assign observer to moisture part of the engine.
    //!
    //! @param observer observer that will be notified about each timestep
    void subscribeMoisture(MultiDomain & domain, Timesteps::TimestepObserver * observer);

    //! \brief Unsubscribe from moisture notifications
    //!
    //! @param observer observer that will be notified about each timestep
    void unsubscribeMoisture(MultiDomain & domain, Timesteps::TimestepObserver * observer);
}   // namespace HygroThermFEM