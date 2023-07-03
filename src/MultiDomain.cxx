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
    MultiDomain::MultiDomain(bool performThermal, bool performMoisture) :
        simulateThermal(performThermal), simulateMoisture(performMoisture)
    {}
}   // namespace HygroThermFEM
