#include "BoundaryCondition2DThermal.hxx"

#include "Common.hxx"
#include "EnclosureRadiation.hxx"
#include "Nodes.hxx"
#include "VectorOperators.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    ConstantConvectionBC::ConstantConvectionBC(Nodes & nodePool,
                                               size_t index1,
                                               size_t index2,
                                               const FixedBCHCCoefficients & fixedBCHCCoefficients,
                                               const bool simulateVaporFluxEnergy) :
        IConvectionBC(nodePool,
                      index1,
                      index2,
                      fixedBCHCCoefficients.AirTemperature,
                      ConvectionModelFactory::createFixedFilmCoefficient(
                        m_Nodes, fixedBCHCCoefficients.ConvectionCoefficient),
                      fixedBCHCCoefficients.AirHumidity,
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// ThermalTARPConvectionBC
    ////////////////////////////////////////////////////////
    ThermalTARPConvectionBC::ThermalTARPConvectionBC(Nodes & nodePool,
                                                     size_t index1,
                                                     size_t index2,
                                                     const TARPCoefficients & varHCCoeff,
                                                     const double surfaceTilt,
                                                     const bool simulateVaporFluxEnergy) :
        IConvectionBC(nodePool,
                      index1,
                      index2,
                      varHCCoeff.AirTemperature,
                      ConvectionModelFactory::createTARPFilmCoefficient(
                        m_Nodes, varHCCoeff.AirTemperature, surfaceTilt),
                      varHCCoeff.AirHumidity,
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// ASHRAEInsideConvectionBC
    ////////////////////////////////////////////////////////
    ASHRAEInsideConvectionBC::ASHRAEInsideConvectionBC(Nodes & nodePool,
                                                       size_t index1,
                                                       size_t index2,
                                                       const ASHRAEInsideCoefficients & coeff,
                                                       double surfaceHeight,
                                                       double surfaceTilt,
                                                       bool simulateVaporFluxEnergy) :
        IConvectionBC(
          nodePool,
          index1,
          index2,
          coeff.AirTemperature,
          ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
            m_Nodes, coeff.AirTemperature, surfaceTilt, surfaceHeight, coeff.AirPressure),
          coeff.AirHumidity,
          simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// ASHRAEOutsideConvectionBC
    ////////////////////////////////////////////////////////
    ASHRAEOutsideConvectionBC::ASHRAEOutsideConvectionBC(Nodes & nodePool,
                                                         size_t index1,
                                                         size_t index2,
                                                         const ASHRAEOutsideCoefficients & coeff,
                                                         bool simulateVaporFluxEnergy) :
        IConvectionBC(
          nodePool,
          index1,
          index2,
          coeff.AirTemperature,
          ConvectionModelFactory::createASHRAEOutsideFilmCoefficient(m_Nodes, coeff.WindSpeed),
          coeff.AirHumidity,
          simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// YazdanianKlemsConvectionBC
    ////////////////////////////////////////////////////////
    YazdanianKlemsConvectionBC::YazdanianKlemsConvectionBC(Nodes & nodePool,
                                                           size_t index1,
                                                           size_t index2,
                                                           const YazdanianKlemsCoefficients & coeff,
                                                           bool simulateVaporFluxEnergy) :
        IConvectionBC(nodePool,
                      index1,
                      index2,
                      coeff.AirTemperature,
                      ConvectionModelFactory::createYazdanianKlemsFilmCoefficient(
                        m_Nodes, coeff.AirTemperature, coeff.WindSpeed, coeff.WindDir),
                      coeff.AirHumidity,
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// KimuraConvectionBC
    ////////////////////////////////////////////////////////
    KimuraConvectionBC::KimuraConvectionBC(Nodes & nodePool,
                                           size_t index1,
                                           size_t index2,
                                           const KimuraCoefficients & coeff,
                                           bool simulateVaporFluxEnergy) :
        IConvectionBC(nodePool,
                      index1,
                      index2,
                      coeff.AirTemperature,
                      ConvectionModelFactory::createKimuraFilmCoefficient(
                        m_Nodes, coeff.WindSpeed, coeff.WindDir),
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    // The huge film coefficient is a NUMERICAL pin, not a physical surface, so the
    // vapor-flux energy term must stay off: that term scales with the same film
    // coefficient (beta = h_c / (rho_air cp_air)), and at 1e18 it books a phantom
    // latent flux of order 1e19 W/m^2 against air at humidity zero -- enough to
    // drag the "fixed" temperature tens of kelvin off its pin and stall the
    // coupled nonlinear iteration against it whenever moisture is simulated.
    TemperatureBC::TemperatureBC(Nodes & nodePool,
                                 const size_t index1,
                                 const size_t index2,
                                 const double t_NodeTemperatures) :
        ConstantConvectionBC(nodePool,
                             index1,
                             index2,
                             {t_NodeTemperatures, hugeFilmCoefficient},
                             /*simulateVaporFluxEnergy=*/false)
    {
        auto & node1 = nodePool.getNode(index1);
        auto & node2 = nodePool.getNode(index2);
        node1.setTemperature(t_NodeTemperatures);
        node2.setTemperature(t_NodeTemperatures);
    }

    TemperatureBC::TemperatureBC(Nodes & nodePool,
                                 const size_t index1,
                                 const size_t index2,
                                 const double t_Temp1,
                                 const double t_Temp2) :
        ConstantConvectionBC(nodePool,
                             index1,
                             index2,
                             {(t_Temp1 + t_Temp2) / 2, hugeFilmCoefficient},
                             /*simulateVaporFluxEnergy=*/false)
    {
        auto & node1 = nodePool.getNode(index1);
        auto & node2 = nodePool.getNode(index2);
        node1.setTemperature(t_Temp1);
        node2.setTemperature(t_Temp2);
    }

    ////////////////////////////////////////////////////////
    /// Flux BC
    ////////////////////////////////////////////////////////

    FluxBC::FluxBC(Nodes & nodePool,
                   const size_t index1,
                   const size_t index2,
                   const double t_Flux) :
        IBCLinear2D(nodePool, index1, index2),
        m_Flux(t_Flux)
    {}

    BCVector FluxBC::R_Vector() const
    {
        return m_PsiVector * m_Flux;
    }

    BCMatrix2D FluxBC::H_Matrix() const
    {
        // Flux boundary conditions do not have H matrix (It is zero)
        return BCMatrix2D{};
    }

    ////////////////////////////////////////////////////////
    /// IRadiationBC - Base class for radiation BCs
    ////////////////////////////////////////////////////////

    IRadiationBC::IRadiationBC(Nodes & nodePool,
                               const size_t index1,
                               const size_t index2,
                               const double radiationTemperature,
                               const bool isLinear) :
        IBCLinear2D(nodePool, index1, index2, isLinear),
        m_RadiationTemperature{radiationTemperature}
    {}

    double IRadiationBC::gaussRadiationCoefficient(const std::size_t integrationPointIndex) const
    {
        return radiationCoefficientAt(
          gaussPointProperty(integrationPointIndex, Variable::temperature));
    }

    BCVector IRadiationBC::R_Vector() const
    {
        const auto coefficientAt = [this](const std::size_t idx)
        {
            return gaussRadiationCoefficient(idx);
        };
        return psiGaussWeighted(coefficientAt) * m_RadiationTemperature;
    }

    BCMatrix2D IRadiationBC::H_Matrix() const
    {
        const auto coefficientAt = [this](const std::size_t idx)
        {
            return gaussRadiationCoefficient(idx);
        };
        return psiPsiGaussWeighted(coefficientAt);
    }

    ////////////////////////////////////////////////////////
    /// BlackBodyRadiationBC
    ////////////////////////////////////////////////////////

    BlackBodyRadiationBC::BlackBodyRadiationBC(Nodes & nodePool,
                                               const size_t index1,
                                               const size_t index2,
                                               const double emissivity,
                                               const double radiationTemperature) :
        IRadiationBC(nodePool, index1, index2, radiationTemperature, false),
        m_Emissivity{emissivity}
    {}

    double BlackBodyRadiationBC::radiationCoefficientAt(const double surfaceTemperature) const
    {
        const double surfaceTemp = celsiusToKelvin(surfaceTemperature);
        const double radiationTemp = celsiusToKelvin(m_RadiationTemperature);
        return (surfaceTemp + radiationTemp)
               * (radiationTemp * radiationTemp + surfaceTemp * surfaceTemp)
               * Constants::STEFANBOLTZMANN * m_Emissivity;
    }

    ////////////////////////////////////////////////////////
    /// LinearizedRadiationBC
    ////////////////////////////////////////////////////////

    LinearizedRadiationBC::LinearizedRadiationBC(Nodes & nodePool,
                                                 const size_t index1,
                                                 const size_t index2,
                                                 const LinearizedRadiationBCCoefficients & linearRadBC) :
        IRadiationBC(nodePool, index1, index2, linearRadBC.RadiationTemperature, true),
        m_RadiationCoefficient(linearRadBC.RadiationCoefficient)
    {}

    double LinearizedRadiationBC::radiationCoefficientAt(const double /*surfaceTemperature*/) const
    {
        return m_RadiationCoefficient;
    }

    ////////////////////////////////////////////////////////
    /// EnclosureRadiationBC
    ////////////////////////////////////////////////////////

    EnclosureRadiationBC::EnclosureRadiationBC(Nodes & nodePool,
                                               const size_t index1,
                                               const size_t index2,
                                               const double emissivity,
                                               EnclosureRadiation & coordinator,
                                               const size_t segmentIndex) :
        IBCLinear2D(nodePool, index1, index2, false),   // nonlinear: depends on temperatures
        m_Emissivity{emissivity},
        m_Coordinator{coordinator},
        m_SegmentIndex{segmentIndex}
    {}

    double
      EnclosureRadiationBC::gaussRadiationCoefficient(const std::size_t integrationPointIndex,
                                                      const double radiantTemperature) const
    {
        const double radiationTemp = celsiusToKelvin(radiantTemperature);
        const double surfaceTemp =
          celsiusToKelvin(gaussPointProperty(integrationPointIndex, Variable::temperature));
        return (surfaceTemp + radiationTemp)
               * (radiationTemp * radiationTemp + surfaceTemp * surfaceTemp)
               * Constants::STEFANBOLTZMANN * m_Emissivity;
    }

    double EnclosureRadiationBC::isothermalCoefficient(const double surfaceTemperature,
                                                       const double radiantTemperature) const
    {
        const double surfaceKelvin = celsiusToKelvin(surfaceTemperature);
        const double radiantKelvin = celsiusToKelvin(radiantTemperature);
        return (surfaceKelvin + radiantKelvin)
               * (radiantKelvin * radiantKelvin + surfaceKelvin * surfaceKelvin)
               * Constants::STEFANBOLTZMANN * m_Emissivity;
    }

    BCVector EnclosureRadiationBC::R_Vector() const
    {
        const double radiantTemperature = m_Coordinator.effectiveRadiantTemperature(m_SegmentIndex);
        if(m_Coordinator.surfaceTemperatureModel() == EnclosureSurfaceTemperature::LocalTemperature)
        {
            const auto coefficientAt = [this, radiantTemperature](const std::size_t idx)
            {
                return gaussRadiationCoefficient(idx, radiantTemperature);
            };
            return psiGaussWeighted(coefficientAt) * radiantTemperature;
        }

        // Conrad-compatible (segment-isothermal): the whole segment radiates uniformly at the
        // fourth-power-mean surface temperature, so the net flux h * (tsurf - Trad) is shared
        // equally by the two nodes. H stays the h * PsiPsi tangent (keeps Newton-Raphson
        // convergent); the difference between the interpolated nodal temperature and the uniform
        // tsurf is lagged into R at the current state, which vanishes into the Conrad fixed point
        // at convergence: flux_i = (H*T - R)_i = psi_i * h * (tsurf - Trad).
        const double surfaceTemperature = m_Coordinator.segmentSurfaceTemperature(m_SegmentIndex);
        const double hCoefficient = isothermalCoefficient(surfaceTemperature, radiantTemperature);
        BCVector result{};
        for(std::size_t row = 0; row < numOfBCNodes; ++row)
        {
            double interpolated{0.0};
            for(std::size_t col = 0; col < numOfBCNodes; ++col)
            {
                interpolated +=
                  m_PsiPsiMatrix(row, col) * m_Nodes[col].property(Variable::temperature);
            }
            result[row] = hCoefficient
                          * (m_PsiVector[row] * (radiantTemperature - surfaceTemperature)
                             + interpolated);
        }
        return result;
    }

    BCMatrix2D EnclosureRadiationBC::H_Matrix() const
    {
        const double radiantTemperature = m_Coordinator.effectiveRadiantTemperature(m_SegmentIndex);
        if(m_Coordinator.surfaceTemperatureModel() == EnclosureSurfaceTemperature::LocalTemperature)
        {
            const auto coefficientAt = [this, radiantTemperature](const std::size_t idx)
            {
                return gaussRadiationCoefficient(idx, radiantTemperature);
            };
            return psiPsiGaussWeighted(coefficientAt);
        }
        const double surfaceTemperature = m_Coordinator.segmentSurfaceTemperature(m_SegmentIndex);
        const double hCoefficient = isothermalCoefficient(surfaceTemperature, radiantTemperature);
        return psiPsiGaussWeighted([hCoefficient](const std::size_t) { return hCoefficient; });
    }
}   // namespace HygroThermFEM