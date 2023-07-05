#pragma once

#include <string>
#include <vector>
#include <algorithm>

#include "Common.hxx"
#include "TimestepObserver.hxx"

namespace Timesteps
{
    //! \brief Observable template for timestep change notification.
    //!
    //! Engine will try to perform simulation within given timestep division. Sometimes that is not
    //! working and the engine will divide timestep into smaller divisions. This class implements
    //! notification for outside world in case someone is interested what the engine is doing.

    class TimestepNotifier
    {
    private:
        std::vector<TimestepObserver *> observers;

    public:
        //! \brief Notification to observers of current statue of the simulation
        //!
        //! \param divisionLevel - Current timestep division level
        //! \param timestepNumber - Current timestep index in current division
        void notify(unsigned divisionLevel, unsigned timestepNumber);

        //! \brief Subscribe for notifications.
        //!
        //! \param observer - Observer that needs to be notified about simulation changes
        void subscribe(TimestepObserver * observer);

        //! \brief Unsubscribe from notifications.
        //!
        //! \param observer - Observer that was notified about simulation changes.
        void unsubscribe(TimestepObserver * observer);
    };

}   // namespace Timesteps