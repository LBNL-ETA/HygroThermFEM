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

}   // namespace HygroThermFEM
