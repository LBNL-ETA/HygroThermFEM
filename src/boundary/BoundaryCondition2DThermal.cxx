#include "BoundaryCondition2DThermal.hxx"

#include "Common.hxx"
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

    std::vector<double> IRadiationBC::R_Vector() const
    {
        return m_PsiVector * radiationCoefficients() * m_RadiationTemperature;
    }

    SquareMatrix IRadiationBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(radiationCoefficients());
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

    std::vector<double> BlackBodyRadiationBC::radiationCoefficients() const
    {
        std::vector<double> result(numOfBCNodes, 0);
        for(std::size_t idx = 0; idx < numOfBCNodes; ++idx)
        {
            const double surfaceTemp = celsiusToKelvin(m_Nodes[idx].property(Variable::temperature));
            const double radiationTemp = celsiusToKelvin(m_RadiationTemperature);
            result[idx] = (surfaceTemp + radiationTemp)
                          * (radiationTemp * radiationTemp + surfaceTemp * surfaceTemp)
                          * Constants::STEFANBOLTZMANN * m_Emissivity;
        }
        return result;
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

    std::vector<double> LinearizedRadiationBC::radiationCoefficients() const
    {
        return std::vector<double>(numOfBCNodes, m_RadiationCoefficient);
    }
}   // namespace HygroThermFEM