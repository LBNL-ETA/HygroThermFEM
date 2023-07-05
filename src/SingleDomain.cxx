#include <cmath>

#include "SingleDomain.hxx"
#include "FEMunique.hxx"
#include "LinearSolver.hxx"
#include "Common.hxx"
#include "FEMMath.hxx"
#include "BoundaryCondition2D.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"
#include "NodePool.hxx"
#include "TimestepData.hxx"
#include "MaterialPool.hxx"

namespace HygroThermFEM
{
    SingleDomain::SingleDomain(DomainType type) :
        domainType(type)
    {}

    std::vector<NodeFlux> SingleDomain::flux() const
    {
        return elements.flux();
    }

    void SingleDomain::setGravityVector(const FenestrationCommon::GravityVector & aGravityVector)
    {
        gravityVector = aGravityVector;
        if(gasCavities.has_value())
        {
            gasCavities->setGravityVector(gravityVector);
        }
    }

    BaseVariable baseVariable(SingleDomain & domain)
    {
        const std::map<DomainType, BaseVariable> baseVariableMap = {
          {DomainType::Thermal, BaseVariable::temperature},
          {DomainType::Moisture, BaseVariable::humidity}};

        return baseVariableMap.at(domain.domainType);
    }
}   // namespace HygroThermFEM
