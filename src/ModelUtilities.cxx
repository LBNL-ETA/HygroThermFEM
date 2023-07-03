#include "ModelUtilities.hxx"
#include "MultiDomain.hxx"

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

    void setGravityVector(MultiDomain & domain,
                          const FenestrationCommon::GravityVector & gravityVector)
    {
        domain.thermalDomain.setGravityVector(gravityVector);
    }

    void clearModel(MultiDomain & domain)
    {
        domain.thermalDomain.clearModel();
        domain.moistureDomain.clearModel();
    }
}   // namespace HygroThermFEM