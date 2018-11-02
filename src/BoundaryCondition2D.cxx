#include <algorithm>
#include <cmath>
#include <functional>

#include "BoundaryCondition2D.hxx"
#include "Common.hxx"
#include "MaterialProperties.hxx"
#include "VectorOperators.hxx"

namespace MoisThermFEM
{
    ////////////////////////////////////////////////////////
    /// ConvectionBC
    ////////////////////////////////////////////////////////

    ConvectionBC::ConvectionBC(const Node2D & t_Node1,
                               const Node2D & t_Node2,
                               const double t_ConvectionCoefficient,
                               const double t_AirTemperature) :
        IBCLinear2D(t_Node1, t_Node2),
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

    TemperatureBC::TemperatureBC(Node2D & t_Node1,
                                 Node2D & t_Node2,
                                 const double t_NodeTemperatures) :
        ConvectionBC(t_Node1, t_Node2, 1e18, t_NodeTemperatures)
    {
        t_Node1.setProperty(Property::temperature, t_NodeTemperatures);
        t_Node2.setProperty(Property::temperature, t_NodeTemperatures);
    }

    TemperatureBC::TemperatureBC(Node2D & t_Node1,
                                 Node2D & t_Node2,
                                 const double t_Temp1,
                                 const double t_Temp2) :
        ConvectionBC(t_Node1, t_Node2, 1e18, (t_Temp1 + t_Temp2) / 2)
    {
        t_Node1.setProperty(Property::temperature, t_Temp1);
        t_Node2.setProperty(Property::temperature, t_Temp2);
    }

    ////////////////////////////////////////////////////////
    /// Flux BC
    ////////////////////////////////////////////////////////

    FluxBC::FluxBC(Node2D & t_Node1, Node2D & t_Node2, const double t_Flux) :
        IBCLinear2D(t_Node1, t_Node2),
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

    BlackBodyRadiationBC::BlackBodyRadiationBC(const Node2D & t_Node1,
                                               const Node2D & t_Node2,
                                               const double t_Emissivity,
                                               const double t_RadiationTemperature) :
        IBCLinear2D(t_Node1, t_Node2, false),
        m_RadiationTemperature{t_RadiationTemperature},
        m_Emissivity{t_Emissivity}
    {}

    std::vector<double> BlackBodyRadiationBC::HRadiative() const
    {
        std::vector<double> result(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            double T = m_Nodes[j].getProperty(Property::temperature);
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

    MoistureBC::MoistureBC(const Node2D & t_Node1,
                           const Node2D & t_Node2,
                           const double t_ConvectiveCoefficient,
                           const Material & t_Material,
                           const double t_AirHumidity,
                           const double t_AirTemperature) :
        IBCLinear2D(t_Node1, t_Node2),
        m_ConvectiveCoefficient(t_ConvectiveCoefficient),
        m_AirHumidity(t_AirHumidity),
        m_AirTemperature(t_AirTemperature),
        m_Material(t_Material)
    {}

    std::vector<double> MoistureBC::R_Vector() const
    {
        // pValue airFill = MaterialProperties::getAirFill( m_Material );
        auto humidityCalculator = SaturationFunction() * m_Material.porosity() * m_AirHumidity;

        const auto humidityByVolume =
          humidityCalculator.value(State(m_AirTemperature, m_AirHumidity, 101325, 0));
        const auto coeff =
          m_ConvectiveCoefficient * humidityByVolume / (Constants::Density_Air * Constants::Cp_Air);
        return m_PsiVector * coeff;
    }

    FenestrationCommon::SquareMatrix MoistureBC::H_Matrix() const
    {
        const auto humidityCoeff = m_Material.porosity() * SaturationFunction();

        const auto humidityByVolume1 = humidityCoeff.value(m_Nodes[0].getState());
        const auto humidityByVolume2 = humidityCoeff.value(m_Nodes[1].getState());

        const auto beta = m_ConvectiveCoefficient / (Constants::Density_Air * Constants::Cp_Air);

        std::vector<double> coeffs{humidityByVolume1 * beta, humidityByVolume2 * beta};

        return m_PsiPsiMatrix.mmultRows(coeffs);
    }
}   // namespace MoisThermFEM