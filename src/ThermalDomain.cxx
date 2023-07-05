#include "ThermalDomain.hxx"

#include "BoundaryCondition2D.hxx"
#include "BoundaryCondition2DThermal.hxx"
#include "FEMunique.hxx"

namespace HygroThermFEM
{
    void ThermalDomain::postProcess(std::vector<double> & solution)
    {
        SingleDomain::postProcess(solution);
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
    }

    ThermalDomain::ThermalDomain(bool automaticUpdatePreviousTimestep) :
        SingleDomain(
          DomainType::Thermal, BaseVariable::temperature, automaticUpdatePreviousTimestep)
    {}
}   // namespace HygroThermFEM
