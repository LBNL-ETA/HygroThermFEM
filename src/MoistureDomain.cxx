#include "MoistureDomain.hxx"

#include "BoundaryCondition2D.hxx"
#include "BoundaryCondition2DMoisture.hxx"

namespace HygroThermFEM
{
    MoistureDomain::MoistureDomain(bool automaticUpdatePreviousTimestep) :
        SingleDomain(DomainType::Moisture, BaseVariable::humidity, automaticUpdatePreviousTimestep)
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
