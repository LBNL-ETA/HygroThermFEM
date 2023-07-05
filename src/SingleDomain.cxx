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
    // Define a map from SingleDomainType to ElementFactory.
    std::map<DomainType, ElementFactory> elementFactoryMap = {
      {DomainType::Thermal,
       [](const size_t index1,
          const size_t index2,
          const size_t index3,
          const size_t index4,
          const std::string & materialName) {
           return std::make_unique<ElementThermalLinear2D>(
             index1, index2, index3, index4, materialName);
       }},
      {DomainType::Moisture,
       [](const size_t index1,
          const size_t index2,
          const size_t index3,
          const size_t index4,
          const std::string & materialName) {
           return std::make_unique<ElementMoistureLinear2D>(
             index1, index2, index3, index4, materialName);
       }}};

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
        domainType(type),
        gasCavities(nullptr)
    {}

    void SingleDomain::createElement(
      size_t index1, size_t index2, size_t index3, size_t index4, const std::string & materialName)
    {
        m_Elements.assignElement(
          elementFactoryMap[domainType](index1, index2, index3, index4, materialName));
    }

    std::vector<NodeFlux> SingleDomain::flux() const
    {
        return m_Elements.flux();
    }

    void SingleDomain::postProcess(std::vector<double> & solution)
    {
        if(gasCavities == nullptr)
        {
            gasCavities = std::make_unique<EquivalentGasCavities>(m_Elements);
            gasCavities->setGravityVector(m_GravityVector);
        }
        gasCavities->update();

        // Domain-specific processing.
        postProcessFuncMap[domainType](solution);
    }

    void SingleDomain::setGravityVector(const FenestrationCommon::GravityVector & gravityVector)
    {
        m_GravityVector = gravityVector;
        if(gasCavities != nullptr)
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
