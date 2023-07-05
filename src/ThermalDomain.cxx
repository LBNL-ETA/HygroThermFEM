#include "ThermalDomain.hxx"

#include "BoundaryCondition2D.hxx"
#include "BoundaryCondition2DThermal.hxx"
#include "FEMunique.hxx"

namespace HygroThermFEM
{
    void ThermalDomain::createElement(const size_t index1,
                                      const size_t index2,
                                      const size_t index3,
                                      const size_t index4,
                                      const std::string & materialName)
    {
        m_Elements.assignElement(
          std::make_unique<ElementThermalLinear2D>(index1, index2, index3, index4, materialName));
    }

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
        SingleDomain(BaseVariable::temperature, automaticUpdatePreviousTimestep)
    {}
}   // namespace HygroThermFEM
