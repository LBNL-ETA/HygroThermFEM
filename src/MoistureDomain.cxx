#include "MoistureDomain.hxx"

#include "BoundaryCondition2D.hxx"
#include "BoundaryCondition2DMoisture.hxx"

namespace HygroThermFEM
{
    void MoistureDomain::createElement(const size_t index1,
                                       const size_t index2,
                                       const size_t index3,
                                       const size_t index4,
                                       const std::string & materialName)
    {
        m_Elements.assignElement(
          std::make_unique<ElementMoistureLinear2D>(index1, index2, index3, index4, materialName));
    }

    MoistureDomain::MoistureDomain(bool automaticUpdatePreviousTimestep) :
        SingleDomain(BaseVariable::humidity, automaticUpdatePreviousTimestep)
    {}

    void MoistureDomain::postProcess(std::vector<double> & solution)
    {
        SingleDomain::postProcess(solution);
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
    }
}
