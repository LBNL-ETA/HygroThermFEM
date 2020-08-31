#include <algorithm>

#include "BoundaryCondition2D.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"
#include "VectorOperators.hxx"
#include "Common.hxx"
#include "SimulationProperties.hxx"
#include "ConvectiveCoefficient.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// IConvectionBC
    ////////////////////////////////////////////////////////

    IConvectionBC::IConvectionBC(size_t index1,
                                 size_t index2,
                                 const double t_AirTemperature,
                                 std::unique_ptr<IConvectiveCoefficient> convModel,
                                 const double t_AirHumidity,
                                 const bool t_SimulateMoisture) :
        IBCLinear2D(index1, index2),
        m_AirTemperature(t_AirTemperature),
        m_ConvectiveCoeffCalc(std::move(convModel)),
        m_AirHumidity(t_AirHumidity),
        m_SimulateVaporFluxEnergy(t_SimulateMoisture)
    {}

    std::vector<double> IConvectionBC::R_Vector() const
    {
        auto rightHandSide = m_ConvectiveCoeffCalc->convectiveCoefficients() * m_AirTemperature;

        const auto excludeHeatOfEvaporation =
          SimulationProperties::Instance().excludeHeatOfEvaporation();
        // Moisture is dumping some energy into domain. However, it is possible that user
        // choose not to simulate moisture in which case energy should not be included in
        // simulation
        if(m_SimulateVaporFluxEnergy && !excludeHeatOfEvaporation)
        {
            // Vapor leaking part is added here
            std::vector<double> vaporLeak(numOfBCNodes, 0);
            for(std::size_t j = 0; j < numOfBCNodes; ++j)
            {
                const double T = m_Nodes[j].property(Variable::temperature);
                const double humidity = m_Nodes[j].property(Variable::humidity);
                vaporLeak[j] =
                  (m_AirHumidity * saturationConcentrationAtTemperature(m_AirTemperature)
                   - humidity * saturationConcentrationAtTemperature(T))
                  * heatOfEvaporation(T);
            }
            const auto vaporFluxEnergy =
              vaporLeak * m_ConvectiveCoeffCalc->waterVaporTransferCoefficient();
            rightHandSide = rightHandSide + vaporFluxEnergy;
        }

        return m_PsiVector * rightHandSide;
    }

    SquareMatrix IConvectionBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(m_ConvectiveCoeffCalc->convectiveCoefficients());
    }

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
    /// TARPConvectionBC
    ////////////////////////////////////////////////////////
    TARPConvectionBC::TARPConvectionBC(size_t index1,
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
      IConvectionBC(index1,
                    index2,
                    coeff.AirTemperature,
                    ConvectionModelFactory::createASHRAEInsideFilmCoefficient(
                      m_Nodes, coeff.AirTemperature, surfaceTilt, surfaceHeight, coeff.AirPressure),
                    coeff.AirHumidity,
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
        IBCLinear2D(index1, index2), m_Flux(t_Flux)
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

    /////////////////////////////////////////////////////
    /// IMoistureBC
    /////////////////////////////////////////////////////

    IMoistureBC::IMoistureBC(const size_t index1,
                             const size_t index2,
                             const std::string & materialName,
                             const double t_AirHumidity,
                             const double t_AirTemperature,
                             std::unique_ptr<IConvectiveCoefficient> model) :
        IBCLinear2D(index1, index2),
        m_AirHumidity(t_AirHumidity),
        m_AirTemperature(t_AirTemperature),
        m_Material(MaterialPool::Instance().material(materialName)),
        m_ConvectiveCoeffCalc(std::move(model))
    {}

    std::vector<double> IMoistureBC::R_Vector() const
    {
        const auto satOutside = saturationConcentrationAtTemperature(m_AirTemperature);
        const auto gconv =
          m_ConvectiveCoeffCalc->waterVaporTransferCoefficient() * satOutside * m_AirHumidity;
        return m_PsiVector * gconv;
    }

    SquareMatrix IMoistureBC::H_Matrix() const
    {
        std::vector<double> concentration(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double T = m_Nodes[j].property(Variable::temperature);
            concentration[j] = saturationConcentrationAtTemperature(T);
        }
        const auto vaporFlux =
          concentration * m_ConvectiveCoeffCalc->waterVaporTransferCoefficient();

        return m_PsiPsiMatrix.mmultRows(vaporFlux);
    }

    /////////////////////////////////////////////////////
    /// MoistureBCTARPHc
    /////////////////////////////////////////////////////
    MoistureBCTARPHc::MoistureBCTARPHc(size_t index1,
                                       size_t index2,
                                       const std::string & materialName,
                                       const TARPCoefficients & varHCCoeff,
                                       double surfaceTilt) :
        IMoistureBC(index1,
                    index2,
                    materialName,
                    varHCCoeff.AirHumidity,
                    varHCCoeff.AirTemperature,
                    ConvectionModelFactory::createTARPFilmCoefficient(
                      m_Nodes, varHCCoeff.AirTemperature, surfaceTilt))
    {}


    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////
    MoistureBCFixedHc::MoistureBCFixedHc(size_t index1,
                                         size_t index2,
                                         const std::string & materialName,
                                         const FixedBCHCCoefficients & fixedBchcCoefficients) :
        IMoistureBC(index1,
                    index2,
                    materialName,
                    fixedBchcCoefficients.AirHumidity,
                    fixedBchcCoefficients.AirTemperature,
                    ConvectionModelFactory::createFixedFilmCoefficient(
                      m_Nodes, fixedBchcCoefficients.ConvectionCoefficient))
    {}

}   // namespace HygroThermFEM
