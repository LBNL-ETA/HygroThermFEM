#include "ThermalDomain.hxx"

#include "BoundaryCondition2D.hxx"
#include "BoundaryCondition2DThermal.hxx"
#include "FEMunique.hxx"

namespace HygroThermFEM
{
    void ThermalDomain::createBC_FixedHc(const size_t index1,
                                         const size_t index2,
                                         const FixedBCHCCoefficients & fixedBCHCCoefficients,
                                         const bool t_CalculateMoisture)
    {
        m_BCs.assignBC(fem::make_unique<ConstantConvectionBC>(
          index1, index2, fixedBCHCCoefficients, t_CalculateMoisture));
    }

    void ThermalDomain::createBC_FixedHc(
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
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_TARPHc(const size_t index1,
                                        const size_t index2,
                                        const TARPCoefficients & varHCCoeff,
                                        const double surfaceTilt,
                                        const bool simulateVaporFluxEnergy)
    {
        m_BCs.assignBC(std::make_unique<ThermalTARPConvectionBC>(
          index1, index2, varHCCoeff, surfaceTilt, simulateVaporFluxEnergy));
    }

    void ThermalDomain::createBC_TARPHc(size_t index1,
                                        size_t index2,
                                        const std::vector<TARPCoefficients> & varHCCoeff,
                                        double surfaceTilt,
                                        bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        for(const auto & coeff : varHCCoeff)
        {
            timestepBCs.push_back(std::make_unique<ThermalTARPConvectionBC>(
              index1, index2, coeff, surfaceTilt, simulateVaporFluxEnergy));
        }
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_ASHRAEInsideHc(size_t index1,
                                                size_t index2,
                                                const ASHRAEInsideCoefficients & coeff,
                                                double surfaceHeight,
                                                double surfaceTilt,
                                                bool simulateVaporFluxEnergy)
    {
        m_BCs.assignBC(std::make_unique<ASHRAEInsideConvectionBC>(
          index1, index2, coeff, surfaceHeight, surfaceTilt, simulateVaporFluxEnergy));
    }

    void ThermalDomain::createBC_ASHRAEInsideHc(size_t index1,
                                                size_t index2,
                                                const std::vector<ASHRAEInsideCoefficients> & coeff,
                                                double surfaceHeight,
                                                double surfaceTilt,
                                                bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<ASHRAEInsideConvectionBC>(
              index1, index2, cf, surfaceHeight, surfaceTilt, simulateVaporFluxEnergy));
        }
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                                 size_t index2,
                                                 const ASHRAEOutsideCoefficients & coeff,
                                                 bool simulateVaporFluxEnergy)
    {
        m_BCs.assignBC(std::make_unique<ASHRAEOutsideConvectionBC>(
          index1, index2, coeff, simulateVaporFluxEnergy));
    }

    void
      ThermalDomain::createBC_ASHRAEOutsideHc(size_t index1,
                                              size_t index2,
                                              const std::vector<ASHRAEOutsideCoefficients> & coeff,
                                              bool simulateVaporFluxEnergy)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        for(const auto & cf : coeff)
        {
            timestepBCs.push_back(std::make_unique<ASHRAEOutsideConvectionBC>(
              index1, index2, cf, simulateVaporFluxEnergy));
        }
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_FixedTemperature(const size_t index1,
                                                  const size_t index2,
                                                  double t_Temp1,
                                                  double t_Temp2)
    {
        m_BCs.assignBC(std::make_unique<TemperatureBC>(index1, index2, t_Temp1, t_Temp2));
    }

    void ThermalDomain::createBC_FixedTemperature(size_t index1,
                                                  size_t index2,
                                                  const std::vector<ConstantBCTemperatures> & temp)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(temp.begin(), temp.end(), [&](const auto & bc) {
            timestepBCs.push_back(
              std::make_unique<TemperatureBC>(index1, index2, bc.Temperature1, bc.Temperature2));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_FixedTemperature(const size_t index1,
                                                  const size_t index2,
                                                  const double t_Temp)
    {
        m_BCs.assignBC(fem::make_unique<TemperatureBC>(index1, index2, t_Temp));
    }

    void ThermalDomain::createBC_FixedTemperature(size_t index1,
                                                  size_t index2,
                                                  std::vector<double> temp)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(temp.begin(), temp.end(), [&](const auto & temperature) {
            timestepBCs.push_back(std::make_unique<TemperatureBC>(index1, index2, temperature));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_FixedFlux(const size_t index1,
                                           const size_t index2,
                                           const double t_Flux)
    {
        m_BCs.assignBC(fem::make_unique<FluxBC>(index1, index2, t_Flux));
    }

    void ThermalDomain::createBC_BlackBodyRadiation(const size_t index1,
                                                    const size_t index2,
                                                    const double t_Emissivity,
                                                    const double t_RadiationTemperature)
    {
        m_BCs.assignBC(fem::make_unique<BlackBodyRadiationBC>(
          index1, index2, t_Emissivity, t_RadiationTemperature));
    }

    void ThermalDomain::createBC_BlackBodyRadiation(
      size_t index1, size_t index2, const std::vector<BlackBodyRadiationBCCoefficients> & radCoeffs)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(radCoeffs.begin(), radCoeffs.end(), [&](const auto & bc) {
            timestepBCs.push_back(std::make_unique<BlackBodyRadiationBC>(
              index1, index2, bc.Emissivity, bc.Temperature));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

    void ThermalDomain::createBC_LinearizedRadiation(
      const size_t index1,
      const size_t index2,
      const LinearizedRadiationBCCoefficients & linearRadBC)
    {
        m_BCs.assignBC(fem::make_unique<LinearizedRadiationBC>(index1, index2, linearRadBC));
    }

    void ThermalDomain::createBC_LinearizedRadiation(
      size_t index1,
      size_t index2,
      const std::vector<LinearizedRadiationBCCoefficients> & linearRadBC)
    {
        std::vector<std::unique_ptr<IBCLinear2D>> timestepBCs;
        std::for_each(linearRadBC.begin(), linearRadBC.end(), [&](const auto & bc) {
            timestepBCs.push_back(std::make_unique<LinearizedRadiationBC>(index1, index2, bc));
        });
        m_BCs.assignTimestepBCs(std::move(timestepBCs));
    }

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
        IDomain::postProcess(solution);
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
        IDomain(BaseVariable::temperature, automaticUpdatePreviousTimestep)
    {}
}   // namespace HygroThermFEM
