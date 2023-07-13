#include "ModelUtilities.hxx"
#include "SingleDomain.hxx"
#include "MultiDomain.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"

namespace HygroThermFEM
{
    // Define a map from SingleDomainType to ElementFactory.
    const std::map<DomainType, ElementFactory> elementFactoryMap = {
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

    void createElement(SingleDomain & domain,
                       size_t index1,
                       size_t index2,
                       size_t index3,
                       size_t index4,
                       const std::string & materialName)
    {
        domain.elements.assignElement(
          elementFactoryMap.at(domain.domainType)(index1, index2, index3, index4, materialName));
    }

    void createElement(MultiDomain & domain,
                       const size_t index1,
                       const size_t index2,
                       const size_t index3,
                       const size_t index4,
                       const std::string & materialName)
    {
        createElement(domain.thermalDomain, index1, index2, index3, index4, materialName);
        createElement(domain.moistureDomain, index1, index2, index3, index4, materialName);
    }

    bool isLinear(const SingleDomain & domain)
    {
        return domain.elements.isLinear() && domain.boundaryConditions.isLinear();
    }

    void clearModel(SingleDomain & domain)
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
        domain.boundaryConditions.clear();
        domain.elements.clearElements();
    }

    void clearModel(MultiDomain & domain)
    {
        clearModel(domain.thermalDomain);
        clearModel(domain.moistureDomain);
    }

    void setGravityVector(MultiDomain & domain,
                          const FenestrationCommon::GravityVector & gravityVector)
    {
        domain.thermalDomain.setGravityVector(gravityVector);
    }
}   // namespace HygroThermFEM