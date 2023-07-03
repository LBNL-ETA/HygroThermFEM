#include <cmath>
#include <limits>
#include <utility>

#include "MultiDomain.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "NodePool.hxx"
#include "SimulationProperties.hxx"
#include "MaterialDataChecker.hxx"

namespace HygroThermFEM
{
    void MultiDomain::createElement(const size_t index1,
                                    const size_t index2,
                                    const size_t index3,
                                    const size_t index4,
                                    const std::string & materialName)
    {
        thermalDomain.createElement(index1, index2, index3, index4, materialName);
        moistureDomain.createElement(index1, index2, index3, index4, materialName);
    }

    std::vector<double> MultiDomain::property(Variable property)
    {
        return NodePool::Instance().properties(property);
    }

    void MultiDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        thermalDomain.setGravityVector(gravityVector);
    }

    void MultiDomain::subscribeThermal(Timesteps::TimestepObserver * observer)
    {
        thermalDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeThermal(Timesteps::TimestepObserver * observer)
    {
        thermalDomain.unsubscribe(observer);
    }

    void MultiDomain::subscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        moistureDomain.subscribe(observer);
    }

    void MultiDomain::unsubscribeMoisture(Timesteps::TimestepObserver * observer)
    {
        moistureDomain.unsubscribe(observer);
    }

    void MultiDomain::clearModel()
    {
        thermalDomain.clearModel();
        moistureDomain.clearModel();
    }

    MultiDomain::MultiDomain(bool performThermal, bool performMoisture) :
        simulateThermal(performThermal), simulateMoisture(performMoisture)
    {}
}   // namespace HygroThermFEM
