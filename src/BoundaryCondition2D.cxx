#include <algorithm>

#include "BoundaryCondition2D.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"
#include "VectorOperators.hxx"
#include "Common.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////
    /// ConvectionBC
    ////////////////////////////////////////////////////////

    ConvectionBC::ConvectionBC(const size_t index1,
                               const size_t index2,
                               const double t_ConvectionCoefficient,
                               const double t_AirTemperature) :
        IBCLinear2D(index1, index2),
        m_ConvectionCoefficient(t_ConvectionCoefficient),
        m_AirTemperature(t_AirTemperature)
    {}

    std::vector<double> ConvectionBC::R_Vector() const
    {
        return m_PsiVector * m_ConvectionCoefficient * m_AirTemperature;
    }

    FenestrationCommon::SquareMatrix ConvectionBC::H_Matrix() const
    {
        return m_PsiPsiMatrix * m_ConvectionCoefficient;
    }

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_NodeTemperatures) :
        ConvectionBC(index1, index2, 1e18, t_NodeTemperatures)
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
        ConvectionBC(index1, index2, 1e18, (t_Temp1 + t_Temp2) / 2)
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
        std::transform(m_PsiVector.begin(),
                       m_PsiVector.end(),
                       result.begin(),
                       std::bind1st(std::multiplies<double>(), m_Flux));
        return result;
    }

    FenestrationCommon::SquareMatrix FluxBC::H_Matrix() const
    {
        // Flux boundary conditions do not have H matrix (It is zero)
        return FenestrationCommon::SquareMatrix(4);
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
            double T = m_Nodes[j].property(Variable::temperature);
            result[j] = (T + m_RadiationTemperature)
                        * (std::pow(T, 2) + std::pow(m_RadiationTemperature, 2))
                        * Constants::STEFANBOLTZMANN * m_Emissivity;
        }
        return result;
    }

    std::vector<double> BlackBodyRadiationBC::R_Vector() const
    {
        return m_PsiVector * HRadiative() * m_RadiationTemperature;
    }

    FenestrationCommon::SquareMatrix BlackBodyRadiationBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(HRadiative());
    }

    /////////////////////////////////////////////////////
    /// MoistureBC
    /////////////////////////////////////////////////////

    MoistureBC::MoistureBC(size_t index1,
                           size_t index2,
                           const std::string & materialName,
                           double t_AirHumidity,
                           double t_AirTemperature) :
        IBCLinear2D(index1, index2),
        m_AirHumidity(t_AirHumidity),
        m_AirTemperature(t_AirTemperature),
        m_Material(MaterialPool::Instance().material(materialName))
    {}

    std::vector<double> MoistureBC::R_Vector() const
    {
        const auto satOutside = boundarySaturationAtTemperature(m_AirTemperature);
        const auto gconv = betaConv() * satOutside * m_AirHumidity;
        return m_PsiVector * gconv;
    }

    FenestrationCommon::SquareMatrix MoistureBC::H_Matrix() const
    {
        std::vector<double> saturation(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            double T = m_Nodes[j].property(Variable::temperature);
            saturation[j] = boundarySaturationAtTemperature(T);
        }
        const auto coeffs = saturation * betaConv();

        return m_PsiPsiMatrix.mmultRows(coeffs);
    }

    std::vector<double> MoistureBC::betaConv() const
    {
        std::vector<double> betaCon(numOfBCNodes, 0);

        auto minimumConvectionCoefficient = 0.0;

        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double humidity = m_Nodes[j].property(Variable::humidity);
            if(humidity <= 1)
            {
                const double T = m_Nodes[j].property(Variable::temperature);
                const auto convectiveCoefficient =
                  std::max(minimumConvectionCoefficient,
                           1.31 * std::pow(std::abs(T - m_AirTemperature), 1.0 / 3.0));
                betaCon[j] = 7e-9 * convectiveCoefficient;
            }
            else
            {
                betaCon[j] = 0;
            }
        }
        return betaCon;
    }
}   // namespace MoisThermFEM
