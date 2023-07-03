#include "SimulationObserver.hxx"
#include "MultiDomain.hxx"

namespace HygroThermFEM
{
    void subscribeThermal(MultiDomain & domain, Timesteps::TimestepObserver * observer)
    {
        domain.thermalDomain.subscribe(observer);
    }

    void unsubscribeThermal(MultiDomain & domain, Timesteps::TimestepObserver * observer)
    {
        domain.thermalDomain.unsubscribe(observer);
    }

    void subscribeMoisture(MultiDomain & domain, Timesteps::TimestepObserver * observer)
    {
        domain.moistureDomain.subscribe(observer);
    }

    void unsubscribeMoisture(MultiDomain & domain, Timesteps::TimestepObserver * observer)
    {
        domain.moistureDomain.unsubscribe(observer);
    }
}   // namespace HygroThermFEM