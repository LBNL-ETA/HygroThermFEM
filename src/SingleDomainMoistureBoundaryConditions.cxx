#include "SingleDomainMoistureBoundaryConditions.hxx"
#include "SingleDomain.hxx"
#include "BoundaryCondition2DMoisture.hxx"

namespace HygroThermFEM::Moisture
{
    void createBC_TARPHc(SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(std::make_unique<MoistureBCTARPHc>(
          index1, index2, Material.name(), varHCCoeff, surfaceTilt));
    }

    void createBC_TARPHc(SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const std::vector<TARPCoefficients> & varCoeff,
                         double surfaceTilt)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        timestepBCs.reserve(varCoeff.size());
        for(const auto & coeff : varCoeff)
        {
            timestepBCs.push_back(std::make_unique<MoistureBCTARPHc>(
              index1, index2, Material.name(), coeff, surfaceTilt));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_ASHRAEInsideHc(SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt)
    {
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(std::make_unique<MoistureBCASHRAEInside>(
          index1, index2, Material.name(), coeff, surfaceHeight, surfaceTilt));
    }

    void createBC_ASHRAEInsideHc(SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const std::vector<ASHRAEInsideCoefficients> & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<MoistureBCASHRAEInside>(
              index1, index2, Material.name(), cf, surfaceHeight, surfaceTilt));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_ASHRAEOutsideHc(SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff)
    {
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(
          std::make_unique<MoistureBCASHRAEOutside>(index1, index2, Material.name(), coeff));
    }

    void createBC_ASHRAEOutsideHc(SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<ASHRAEOutsideCoefficients> & coeff)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(
              std::make_unique<MoistureBCASHRAEOutside>(index1, index2, Material.name(), cf));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_YazdanianKlemsHc(SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff)
    {
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(
          std::make_unique<MoistureYazdanianKlemsBC>(index1, index2, Material.name(), coeff));
    }

    void createBC_YazdanianKlemsHc(SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<YazdanianKlemsCoefficients> & coeff)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(
              std::make_unique<MoistureYazdanianKlemsBC>(index1, index2, Material.name(), cf));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_KimuraHc(SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff)
    {
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(
          std::make_unique<MoistureKimuraBC>(index1, index2, Material.name(), coeff));
    }

    void createBC_KimuraHc(SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const std::vector<KimuraCoefficients> & coeff)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        const auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(
              std::make_unique<MoistureKimuraBC>(index1, index2, Material.name(), cf));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_FixedHc(SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const FixedBCHCCoefficients & fixedBchcCoefficients)
    {
        /// Need to pull material for current moisture boundary condition
        auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(std::make_unique<MoistureBCFixedHc>(
          index1, index2, Material.name(), fixedBchcCoefficients));
    }

    void createBC_FixedHc(SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const std::vector<FixedBCHCCoefficients> & fixedBchcCoefficients)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        std::for_each(
          fixedBchcCoefficients.begin(), fixedBchcCoefficients.end(), [&](const auto & bc) {
              timestepBCs.push_back(
                std::make_unique<MoistureBCFixedHc>(index1, index2, Material.name(), bc));
          });
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_FixedHumidity(SingleDomain & domain,
                                size_t index1,
                                size_t index2,
                                const TemperatureAndHumidity & values)
    {
        auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        domain.m_BCs.assignBC(
          std::make_unique<MoistureBCFixedHumidity>(index1, index2, Material.name(), values));
    }

    void createBC_FixedHumidity(SingleDomain & domain,
                                size_t index1,
                                size_t index2,
                                const std::vector<TemperatureAndHumidity> & values)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs{values.size()};
        auto & Material = domain.m_Elements.findElement(index1, index2)->getMaterial();
        for(size_t i = 0u; i < values.size(); ++i)
        {
            timestepBCs[i] =
              std::make_unique<MoistureBCFixedHumidity>(index1, index2, Material.name(), values[i]);
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }
}   // namespace HygroThermFEM::Moisture