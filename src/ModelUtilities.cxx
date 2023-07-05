#include "ModelUtilities.hxx"
#include "SingleDomain.hxx"
#include "MultiDomain.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"

namespace HygroThermFEM
{
    void createElement(MultiDomain & domain,
                       const size_t index1,
                       const size_t index2,
                       const size_t index3,
                       const size_t index4,
                       const std::string & materialName)
    {
        domain.thermalDomain.createElement(index1, index2, index3, index4, materialName);
        domain.moistureDomain.createElement(index1, index2, index3, index4, materialName);
    }

    bool isLinear(SingleDomain & domain)
    {
        return domain.m_Elements.isLinear() && domain.m_BCs.isLinear();
    }

    void HygroThermFEM::clearModel(SingleDomain & domain)
    {
        NodePool::Instance().clear();
        MaterialPool::Instance().clear();
        domain.m_BCs.clear();
        domain.m_Elements.clearElements();
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