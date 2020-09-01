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

    void MoistureDomain::createBC_TARPHc(const size_t index1,
                                         const size_t index2,
                                         const TARPCoefficients & varHCCoeff,
                                         double surfaceTilt)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(std::make_unique<MoistureBCTARPHc>(
          index1, index2, Material.name(), varHCCoeff, surfaceTilt));
    }

    void MoistureDomain::createBC_TARPHc(size_t index1,
                                         size_t index2,
                                         const std::vector<TARPCoefficients> & varCoeff,
                                         double surfaceTilt)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        for(const auto & coeff : varCoeff)
        {
            timestepBCs.push_back(std::make_unique<MoistureBCTARPHc>(
              index1, index2, Material.name(), coeff, surfaceTilt));
        }
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void MoistureDomain::createBC_ASHRAEInsideHc(size_t index1,
        size_t index2,
        const ASHRAEInsideCoefficients & coeff,
        double surfaceHeight,
        double surfaceTilt)
    {
        const auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(std::make_unique<MoistureBCASHRAEInside>(
          index1, index2, Material.name(), coeff, surfaceHeight, surfaceTilt));
    }

    void MoistureDomain::createBC_ASHRAEInsideHc(size_t index1,
        size_t index2,
        const std::vector<ASHRAEInsideCoefficients> & coeff,
        double surfaceHeight,
        double surfaceTilt)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<MoistureBCASHRAEInside>(
              index1, index2, Material.name(), cf, surfaceHeight, surfaceTilt));
        }
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void MoistureDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                                  size_t index2,
                                                  const ASHRAEOutsideCoefficients & coeff)
    {
        const auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(std::make_unique<MoistureBCASHRAEOutside>(
          index1, index2, Material.name(), coeff));
    }

    void MoistureDomain::createBC_ASHRAEOutsideHc(
      size_t index1, size_t index2, const std::vector<ASHRAEOutsideCoefficients> & coeff)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<MoistureBCASHRAEOutside>(
              index1, index2, Material.name(), cf));
        }
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void MoistureDomain::createBC_FixedHc(const size_t index1,
                                          const size_t index2,
                                          const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        m_BCs.assignBC(std::make_unique<MoistureBCFixedHc>(
          index1, index2, Material.name(), fixedBchcCoefficients));
    }

    void MoistureDomain::createBC_FixedHc(
      size_t index1,
      size_t index2,
      const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        auto & Material = m_Elements.findElement(index1, index2)->getMaterial();
        std::for_each(
          fixedBchcCoefficients.begin(), fixedBchcCoefficients.end(), [&](const auto & bc) {
              timestepBCs.push_back(
                std::make_unique<MoistureBCFixedHc>(index1, index2, Material.name(), bc));
          });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    MoistureDomain::MoistureDomain(bool automaticUpdatePreviousTimestep) :
        IDomain(BaseVariable::humidity, automaticUpdatePreviousTimestep)
    {}

    void MoistureDomain::postProcess(std::vector<double> & solution)
    {
        IDomain::postProcess(solution);
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
