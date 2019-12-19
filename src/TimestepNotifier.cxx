#include "TimestepNotifier.hxx"

namespace Timesteps
{
    void TimestepNotifier::notify(const unsigned divisionLevel, const unsigned timestepNumber)
    {
        for(auto observer : observers)
        {
            observer->levelChanged(divisionLevel, timestepNumber);
        }
    }

    void TimestepNotifier::subscribe(TimestepObserver * observer)
    {
        observers.push_back(observer);
    }

    void TimestepNotifier::unsubscribe(TimestepObserver * observer)
    {
        observers.erase(std::remove(std::begin(observers), std::end(observers), observer),
                        std::end(observers));
    }
}   // namespace Timesteps
