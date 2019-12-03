#pragma once

#include <string>
#include <vector>
#include <algorithm>

#include "Common.hxx"

template<typename>
class TimestepObserver;

//! \brief Observable template for timestep change notification.
//!
//! Engine will try to perform simulation within given timestep division. Sometimes that is not
//! working and the engine will divide timestep into smaller divisions. This class implements
//! notification for outside world in case someone is interested what the engine is doing.
//! \tparam T Observer that needs to implement function which will be used to perform notification
//! actions.

template<typename T>
class TimestepNotifier
{
private:
    std::vector<TimestepObserver<T> *> observers;

protected:
    void notify(const Timestep::Level& timestepLevel, const unsigned timestepNumber)
    {
        for(auto observer : observers)
        {
            observer->levelChanged(timestepLevel, timestepNumber);
        }
    }

public:
    //! \brief Subscribe for notificiations.
    void subscribe(TimestepObserver<T> * observer)
    {
        observers.push_back(observer);
    }

    //! \brief Call itself to unsubscribe from notifications.
    void unsubscribe(TimestepObserver<T> * observer)
    {
        observers.erase(std::remove(begin(observers), end(observers), observer), end(observers));
    }
};
