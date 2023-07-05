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
    // Define a map from SingleDomainType to PostProcessFunc.
    std::map<DomainType, PostProcessFunc> postProcessFuncMap = {
      {DomainType::Thermal,
       [](std::vector<double> & solution) {
           for(auto & val : solution)
           {
               if(val < Constants::ABSOLUTEZERO)
               {
                   val = Constants::ABSOLUTEZERO + 1e-6;
               }
               if(val > 1000)
               {
                   val = 1000;
               }
           }
       }},
      {DomainType::Moisture, [](std::vector<double> & solution) {
           for(auto & val : solution)
           {
               if(val > 1)
               {
                   val = 1;
               }
               if(val < 0)
               {
                   val = 0;
               }
           }
       }}};

    SingleDomain::SingleDomain(DomainType type) :
        domainType(type)
    {}

    std::vector<NodeFlux> SingleDomain::flux() const
    {
        return m_Elements.flux();
    }

    void SingleDomain::postProcess(std::vector<double> & solution)
    {
        if(!gasCavities.has_value())
        {
            gasCavities.emplace(m_Elements);
            gasCavities->setGravityVector(m_GravityVector);
        }
        gasCavities->update();

        // Domain-specific processing.
        postProcessFuncMap[domainType](solution);
    }

    void SingleDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        m_GravityVector = gravityVector;
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
