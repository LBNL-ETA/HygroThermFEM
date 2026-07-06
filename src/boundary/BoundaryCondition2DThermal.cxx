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

    TemperatureBC::TemperatureBC(Nodes & nodePool,
                                 const size_t index1,
                                 const size_t index2,
                                 const double t_NodeTemperatures) :
        ConstantConvectionBC(nodePool, index1, index2, {t_NodeTemperatures, hugeFilmCoefficient})
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
        ConstantConvectionBC(nodePool, index1, index2, {(t_Temp1 + t_Temp2) / 2, hugeFilmCoefficient})
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

    std::vector<double> FluxBC::R_Vector() const
    {
        std::vector<double> result(m_PsiVector.size(), 0);
        std::transform(m_PsiVector.begin(), m_PsiVector.end(), result.begin(), [&](auto && data) {
            return data * m_Flux;
        });
        return result;
    }

    SquareMatrix FluxBC::H_Matrix() const
    {
        // Flux boundary conditions do not have H matrix (It is zero)
        return SquareMatrix(4);
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

    std::vector<double> IRadiationBC::gaussRadiationCoefficients() const
    {
        std::vector<double> result(numOfIntegrationPoints(), 0.0);
        for(std::size_t idx = 0; idx < result.size(); ++idx)
        {
            result[idx] = radiationCoefficientAt(gaussPointProperty(idx, Variable::temperature));
        }
        return result;
    }

    std::vector<double> IRadiationBC::R_Vector() const
    {
        return psiGaussWeighted(gaussRadiationCoefficients()) * m_RadiationTemperature;
    }

    SquareMatrix IRadiationBC::H_Matrix() const
    {
        return psiPsiGaussWeighted(gaussRadiationCoefficients());
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

    std::vector<double>
      EnclosureRadiationBC::gaussRadiationCoefficients(const double radiantTemperature) const
    {
        std::vector<double> result(numOfIntegrationPoints(), 0.0);
        const double radiationTemp = celsiusToKelvin(radiantTemperature);
        for(std::size_t idx = 0; idx < result.size(); ++idx)
        {
            const double surfaceTemp =
              celsiusToKelvin(gaussPointProperty(idx, Variable::temperature));
            result[idx] = (surfaceTemp + radiationTemp)
                          * (radiationTemp * radiationTemp + surfaceTemp * surfaceTemp)
                          * Constants::STEFANBOLTZMANN * m_Emissivity;
        }
        return result;
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

    std::vector<double> EnclosureRadiationBC::R_Vector() const
    {
        const double radiantTemperature = m_Coordinator.effectiveRadiantTemperature(m_SegmentIndex);
        if(m_Coordinator.surfaceTemperatureModel() == EnclosureSurfaceTemperature::LocalTemperature)
        {
            return psiGaussWeighted(gaussRadiationCoefficients(radiantTemperature))
                   * radiantTemperature;
        }

        // Conrad-compatible (segment-isothermal): the whole segment radiates uniformly at the
        // fourth-power-mean surface temperature, so the net flux h * (tsurf - Trad) is shared
        // equally by the two nodes. H stays the h * PsiPsi tangent (keeps Newton-Raphson
        // convergent); the difference between the interpolated nodal temperature and the uniform
        // tsurf is lagged into R at the current state, which vanishes into the Conrad fixed point
        // at convergence: flux_i = (H*T - R)_i = psi_i * h * (tsurf - Trad).
        const double surfaceTemperature = m_Coordinator.segmentSurfaceTemperature(m_SegmentIndex);
        const double hCoefficient = isothermalCoefficient(surfaceTemperature, radiantTemperature);
        std::vector<double> result(numOfBCNodes, 0.0);
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

    SquareMatrix EnclosureRadiationBC::H_Matrix() const
    {
        const double radiantTemperature = m_Coordinator.effectiveRadiantTemperature(m_SegmentIndex);
        if(m_Coordinator.surfaceTemperatureModel() == EnclosureSurfaceTemperature::LocalTemperature)
        {
            return psiPsiGaussWeighted(gaussRadiationCoefficients(radiantTemperature));
        }
        const double surfaceTemperature = m_Coordinator.segmentSurfaceTemperature(m_SegmentIndex);
        const std::vector<double> gaussCoefficients(
          numOfIntegrationPoints(), isothermalCoefficient(surfaceTemperature, radiantTemperature));
        return psiPsiGaussWeighted(gaussCoefficients);
    }
}   // namespace HygroThermFEM