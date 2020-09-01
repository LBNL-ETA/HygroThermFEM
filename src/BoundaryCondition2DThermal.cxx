#include "BoundaryCondition2DThermal.hxx"

#include "Common.hxx"
#include "NodePool.hxx"
#include "VectorOperators.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    ConstantConvectionBC::ConstantConvectionBC(size_t index1,
                                               size_t index2,
                                               const FixedBCHCCoefficients & fixedBCHCCoefficients,
                                               const bool simulateVaporFluxEnergy) :
        IConvectionBC(index1,
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
    ThermalTARPConvectionBC::ThermalTARPConvectionBC(size_t index1,
                                                     size_t index2,
                                                     const TARPCoefficients & varHCCoeff,
                                                     const double surfaceTilt,
                                                     const bool simulateVaporFluxEnergy) :
        IConvectionBC(index1,
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
    ASHRAEInsideConvectionBC::ASHRAEInsideConvectionBC(size_t index1,
                                                       size_t index2,
                                                       const ASHRAEInsideCoefficients & coeff,
                                                       double surfaceHeight,
                                                       double surfaceTilt,
                                                       bool simulateVaporFluxEnergy) :
        IConvectionBC(
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
    ASHRAEOutsideConvectionBC::ASHRAEOutsideConvectionBC(size_t index1,
                                                         size_t index2,
                                                         const ASHRAEOutsideCoefficients & coeff,
                                                         bool simulateVaporFluxEnergy) :
        IConvectionBC(
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
    YazdanianKlemsConvectionBC::YazdanianKlemsConvectionBC(size_t index1,
                                                           size_t index2,
                                                           const YazdanianKlemsCoefficients & coeff,
                                                           bool simulateVaporFluxEnergy) :
        IConvectionBC(index1,
                      index2,
                      coeff.AirTemperature,
                      ConvectionModelFactory::createYazdanianKlemsFilmCoefficient(
                        m_Nodes, coeff.AirTemperature, coeff.WindSpeed, coeff.WindDir),
                      coeff.AirHumidity,
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// YazdanianKlemsConvectionBC
    ////////////////////////////////////////////////////////
    KimuraConvectionBC::KimuraConvectionBC(size_t index1,
                                           size_t index2,
                                           const KimuraCoefficients & coeff,
                                           bool simulateVaporFluxEnergy) :
        IConvectionBC(index1,
                      index2,
                      coeff.AirTemperature,
                      ConvectionModelFactory::createKimuraFilmCoefficient(
                        m_Nodes, coeff.WindSpeed, coeff.WindDir),
                      simulateVaporFluxEnergy)
    {}

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_NodeTemperatures) :
        ConstantConvectionBC(index1, index2, {t_NodeTemperatures, 1e18})
    {
        auto & node1 = NodePool::Instance().getNode(index1);
        auto & node2 = NodePool::Instance().getNode(index2);
        node1.setStateProperty(BaseVariable::temperature, t_NodeTemperatures);
        node2.setStateProperty(BaseVariable::temperature, t_NodeTemperatures);
    }

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_Temp1,
                                 const double t_Temp2) :
        ConstantConvectionBC(index1, index2, {(t_Temp1 + t_Temp2) / 2, 1e18})
    {
        auto & node1 = NodePool::Instance().getNode(index1);
        auto & node2 = NodePool::Instance().getNode(index2);
        node1.setStateProperty(BaseVariable::temperature, t_Temp1);
        node2.setStateProperty(BaseVariable::temperature, t_Temp2);
    }

    ////////////////////////////////////////////////////////
    /// Flux BC
    ////////////////////////////////////////////////////////

    FluxBC::FluxBC(const size_t index1, const size_t index2, const double t_Flux) :
        IBCLinear2D(index1, index2),
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
    /// BlackBodyRadiationBC
    ////////////////////////////////////////////////////////

    BlackBodyRadiationBC::BlackBodyRadiationBC(const size_t index1,
                                               const size_t index2,
                                               const double t_Emissivity,
                                               const double t_RadiationTemperature) :
        IBCLinear2D(index1, index2, false),
        m_RadiationTemperature{t_RadiationTemperature},
        m_Emissivity{t_Emissivity}
    {}

    std::vector<double> BlackBodyRadiationBC::HRadiative() const
    {
        std::vector<double> result(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double T = celsiusToKelvin(m_Nodes[j].property(Variable::temperature));
            const double Trad = celsiusToKelvin(m_RadiationTemperature);
            result[j] =
              (T + Trad) * (Trad * Trad + T * T) * Constants::STEFANBOLTZMANN * m_Emissivity;
        }
        return result;
    }

    std::vector<double> BlackBodyRadiationBC::R_Vector() const
    {
        return m_PsiVector * HRadiative() * m_RadiationTemperature;
    }

    SquareMatrix BlackBodyRadiationBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(HRadiative());
    }

    ////////////////////////////////////////////////////////
    /// SimplifiedRadiationBC
    ////////////////////////////////////////////////////////

    LinearizedRadiationBC::LinearizedRadiationBC(
      size_t index1, size_t index2, const LinearizedRadiationBCCoefficients & linearRadBC) :
        IBCLinear2D(index1, index2),
        m_RadiationCoefficient(linearRadBC.RadiationCoefficient),
        m_RadiationTemperature(linearRadBC.RadiationTemperature)
    {}

    std::vector<double> LinearizedRadiationBC::R_Vector() const
    {
        return m_PsiVector * m_RadiationCoefficient * m_RadiationTemperature;
    }

    SquareMatrix LinearizedRadiationBC::H_Matrix() const
    {
        return m_PsiPsiMatrix * m_RadiationCoefficient;
    }
}   // namespace HygroThermFEM