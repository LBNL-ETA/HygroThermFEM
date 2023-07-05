#include "SingleDomainThermalBoundaryConditions.hxx"
#include "SingleDomain.hxx"
#include "BoundaryCondition2DThermal.hxx"

namespace HygroThermFEM::Thermal
{
    void createBC_FixedHc(HygroThermFEM::SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const FixedBCHCCoefficients & fixedBCHCCoefficients,
                          bool t_CalculateMoisture)
    {
        domain.m_BCs.assignBC(std::make_unique<ConstantConvectionBC>(
          index1, index2, fixedBCHCCoefficients, t_CalculateMoisture));
    }

    void createBC_FixedHc(HygroThermFEM::SingleDomain & domain,
                          size_t index1,
                          size_t index2,
                          const std::vector<FixedBCHCCoefficients> & fixedBCHCCoefficients,
                          bool calculateMoisture)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(
          fixedBCHCCoefficients.begin(), fixedBCHCCoefficients.end(), [&](const auto & bc) {
              timestepBCs.push_back(
                std::make_unique<ConstantConvectionBC>(index1, index2, bc, calculateMoisture));
          });
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_TARPHc(HygroThermFEM::SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const TARPCoefficients & varHCCoeff,
                         double surfaceTilt,
                         bool simulateVaporFluxEnergy)
    {
        domain.m_BCs.assignBC(std::make_unique<ThermalTARPConvectionBC>(
          index1, index2, varHCCoeff, surfaceTilt, simulateVaporFluxEnergy));
    }

    void createBC_TARPHc(HygroThermFEM::SingleDomain & domain,
                         size_t index1,
                         size_t index2,
                         const std::vector<TARPCoefficients> & varHCCoeff,
                         double surfaceTilt,
                         bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        timestepBCs.reserve(varHCCoeff.size());
        for(const auto & coeff : varHCCoeff)
        {
            timestepBCs.push_back(std::make_unique<ThermalTARPConvectionBC>(
              index1, index2, coeff, surfaceTilt, simulateVaporFluxEnergy));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_ASHRAEInsideHc(HygroThermFEM::SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const ASHRAEInsideCoefficients & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt,
                                 bool simulateVaporFluxEnergy)
    {
        domain.m_BCs.assignBC(std::make_unique<ASHRAEInsideConvectionBC>(
          index1, index2, coeff, surfaceHeight, surfaceTilt, simulateVaporFluxEnergy));
    }

    void createBC_ASHRAEInsideHc(HygroThermFEM::SingleDomain & domain,
                                 size_t index1,
                                 size_t index2,
                                 const std::vector<ASHRAEInsideCoefficients> & coeff,
                                 double surfaceHeight,
                                 double surfaceTilt,
                                 bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<ASHRAEInsideConvectionBC>(
              index1, index2, cf, surfaceHeight, surfaceTilt, simulateVaporFluxEnergy));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_ASHRAEOutsideHc(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const ASHRAEOutsideCoefficients & coeff,
                                  bool simulateVaporFluxEnergy)
    {
        domain.m_BCs.assignBC(std::make_unique<ASHRAEOutsideConvectionBC>(
          index1, index2, coeff, simulateVaporFluxEnergy));
    }

    void createBC_ASHRAEOutsideHc(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<ASHRAEOutsideCoefficients> & coeff,
                                  bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<ASHRAEOutsideConvectionBC>(
              index1, index2, cf, simulateVaporFluxEnergy));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_YazdanianKlemsHc(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const YazdanianKlemsCoefficients & coeff,
                                   bool simulateVaporFluxEnergy)
    {
        domain.m_BCs.assignBC(std::make_unique<YazdanianKlemsConvectionBC>(
          index1, index2, coeff, simulateVaporFluxEnergy));
    }

    void createBC_YazdanianKlemsHc(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<YazdanianKlemsCoefficients> & coeff,
                                   bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<YazdanianKlemsConvectionBC>(
              index1, index2, cf, simulateVaporFluxEnergy));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_KimuraHc(HygroThermFEM::SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const KimuraCoefficients & coeff,
                           bool simulateVaporFluxEnergy)
    {
        domain.m_BCs.assignBC(
          std::make_unique<KimuraConvectionBC>(index1, index2, coeff, simulateVaporFluxEnergy));
    }

    void createBC_KimuraHc(HygroThermFEM::SingleDomain & domain,
                           size_t index1,
                           size_t index2,
                           const std::vector<KimuraCoefficients> & coeff,
                           bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        timestepBCs.reserve(coeff.size());
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(
              std::make_unique<KimuraConvectionBC>(index1, index2, cf, simulateVaporFluxEnergy));
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   double t_Temp1,
                                   double t_Temp2)
    {
        domain.m_BCs.assignBC(std::make_unique<TemperatureBC>(index1, index2, t_Temp1, t_Temp2));
    }

    void createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   const std::vector<ConstantBCTemperatures> & temp)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(temp.begin(), temp.end(), [&](const auto & bc) {
            timestepBCs.push_back(
              std::make_unique<TemperatureBC>(index1, index2, bc.Temperature1, bc.Temperature2));
        });
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   double t_Temp)
    {
        domain.m_BCs.assignBC(std::make_unique<TemperatureBC>(index1, index2, t_Temp));
    }

    void createBC_FixedTemperature(HygroThermFEM::SingleDomain & domain,
                                   size_t index1,
                                   size_t index2,
                                   std::vector<double> t_Temp)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs{t_Temp.size()};
        for(size_t i = 0u; i < t_Temp.size(); ++i)
        {
            timestepBCs[i] = std::make_unique<TemperatureBC>(index1, index2, t_Temp[i]);
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_FixedFlux(HygroThermFEM::SingleDomain & domain,
                            size_t index1,
                            size_t index2,
                            double t_Flux)
    {
        domain.m_BCs.assignBC(std::make_unique<FluxBC>(index1, index2, t_Flux));
    }

    void createBC_FixedFlux(HygroThermFEM::SingleDomain & domain,
                            size_t index1,
                            size_t index2,
                            const std::vector<double> & t_Flux)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs{t_Flux.size()};
        for(size_t i = 0u; i < t_Flux.size(); ++i)
        {
            timestepBCs[i] = std::make_unique<TemperatureBC>(index1, index2, t_Flux[i]);
        }
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_BlackBodyRadiation(HygroThermFEM::SingleDomain & domain,
                                     size_t index1,
                                     size_t index2,
                                     double t_Emissivity,
                                     double t_RadiationTemperature)
    {
        domain.m_BCs.assignBC(std::make_unique<BlackBodyRadiationBC>(
          index1, index2, t_Emissivity, t_RadiationTemperature));
    }

    void
      createBC_BlackBodyRadiation(HygroThermFEM::SingleDomain & domain,
                                  size_t index1,
                                  size_t index2,
                                  const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(radCoeffs.begin(), radCoeffs.end(), [&](const auto & bc) {
            timestepBCs.push_back(std::make_unique<BlackBodyRadiationBC>(
              index1, index2, bc.Emissivity, bc.Temperature));
        });
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void createBC_LinearizedRadiation(HygroThermFEM::SingleDomain & domain,
                                      size_t index1,
                                      size_t index2,
                                      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        domain.m_BCs.assignBC(std::make_unique<LinearizedRadiationBC>(index1, index2, linearRadBC));
    }

    void createBC_LinearizedRadiation(
      HygroThermFEM::SingleDomain & domain,
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(linearRadBC.begin(), linearRadBC.end(), [&](const auto & bc) {
            timestepBCs.push_back(std::make_unique<LinearizedRadiationBC>(index1, index2, bc));
        });
        domain.m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }
}   // namespace HygroThermFEM::Thermal
