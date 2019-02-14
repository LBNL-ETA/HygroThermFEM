#include <algorithm>

#include "BoundaryCondition2D.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"
#include "VectorOperators.hxx"
#include "Common.hxx"

namespace HygroThermFEM
{
    ////////////////////////////////////////////////////////
    /// VariableConvectionCoefficient
    ////////////////////////////////////////////////////////
    VariableConvectionCoefficient::VariableConvectionCoefficient() :
        IConvectiveCoefficient()
    {}

    std::vector<double> VariableConvectionCoefficient::value(const INodes & nodes, const double ambientTemperature) const
    {
        auto minimumConvectionCoefficient = 3.0;
        std::vector<double> result;
        for(const auto & temperature : nodes.properties(Variable::temperature))
        {
            result.push_back(
              std::max(minimumConvectionCoefficient,
                       1.31 * std::pow(std::abs(temperature - ambientTemperature), 1.0 / 3.0)));
        }
        return result;
    }

    ////////////////////////////////////////////////////////
    /// FixedConvectionCoefficient
    ////////////////////////////////////////////////////////
    FixedConvectionCoefficient::FixedConvectionCoefficient() :
        IConvectiveCoefficient()
    {}

    std::vector<double> FixedConvectionCoefficient::value(const INodes & nodes, const double convectiveCoefficient) const
    {
        std::vector<double> result(nodes.size(), convectiveCoefficient);

        return result;
    }

    ////////////////////////////////////////////////////////
    /// ConvectionModelFactory
    ////////////////////////////////////////////////////////
    std::unique_ptr<IConvectiveCoefficient>
    ConvectionModelFactory::create(const ConvectionModel model)
    {
        switch(model)
        {
            case ConvectionModel::Fixed:
                return std::unique_ptr<FixedConvectionCoefficient>(new FixedConvectionCoefficient());
            case ConvectionModel::Variable:
                return std::unique_ptr<VariableConvectionCoefficient>(new VariableConvectionCoefficient());
        }
        return nullptr;
    }

    ////////////////////////////////////////////////////////
    /// IConvectionBC
    ////////////////////////////////////////////////////////

    IConvectionBC::IConvectionBC(size_t index1, size_t index2, double t_AirTemperature, ConvectionModel model)
        :
        IBCLinear2D(index1, index2),
        m_AirTemperature(t_AirTemperature),
        m_ConvectiveCoeffCalc(ConvectionModelFactory::create(model))
    {}

    std::vector<double> IConvectionBC::R_Vector() const
    {
        return m_PsiVector * convectionCoefficients() * m_AirTemperature;
    }

    FenestrationCommon::SquareMatrix IConvectionBC::H_Matrix() const
    {
        return m_PsiPsiMatrix.mmultRows(convectionCoefficients());
    }

    ////////////////////////////////////////////////////////
    /// ConstantConvectionBC
    ////////////////////////////////////////////////////////
    ConstantConvectionBC::ConstantConvectionBC(size_t index1, size_t index2,
                                               double t_AirTemperature,
                                               const double m_ConvectionCoefficient)
        : IConvectionBC(index1, index2, t_AirTemperature, ConvectionModel::Fixed),
          m_ConvectionCoefficient(m_ConvectionCoefficient)
    {}

    std::vector<double> ConstantConvectionBC::convectionCoefficients() const
    {
        return m_ConvectiveCoeffCalc->value(m_Nodes, m_ConvectionCoefficient);
    }

    ////////////////////////////////////////////////////////
    /// VariableConvectionBC
    ////////////////////////////////////////////////////////
    VariableConvectionBC::VariableConvectionBC(size_t index1, size_t index2,
                                               double t_AirTemperature)
        : IConvectionBC(index1, index2, t_AirTemperature, ConvectionModel::Variable)
    {}

    std::vector<double> VariableConvectionBC::convectionCoefficients() const
    {
        return m_ConvectiveCoeffCalc->value(m_Nodes, m_AirTemperature);
    }

    ////////////////////////////////////////////////////////
    /// TemperatureBC
    ////////////////////////////////////////////////////////

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_NodeTemperatures) : ConstantConvectionBC(index1, index2,
                                                                                  t_NodeTemperatures, 1e18)
    {
        auto & node1 = NodePool::Instance().getNode(index1);
        auto & node2 = NodePool::Instance().getNode(index2);
        node1.setStateProperty(BaseVariable::temperature, t_NodeTemperatures);
        node2.setStateProperty(BaseVariable::temperature, t_NodeTemperatures);
    }

    TemperatureBC::TemperatureBC(const size_t index1,
                                 const size_t index2,
                                 const double t_Temp1,
                                 const double t_Temp2) : ConstantConvectionBC(index1, index2,
                                                                       (t_Temp1 + t_Temp2) / 2, 1e18)
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
    /// IMoistureBC
    /////////////////////////////////////////////////////

    IMoistureBC::IMoistureBC(size_t index1,
                             size_t index2,
                             const std::string & materialName,
                             double t_AirHumidity,
                             double t_AirTemperature,
                             ConvectionModel model) :
        IBCLinear2D(index1, index2),
        m_AirHumidity(t_AirHumidity),
        m_AirTemperature(t_AirTemperature),
        m_Material(MaterialPool::Instance().material(materialName)),
        m_ConvectiveCoeffCalc(ConvectionModelFactory::create(model))
    {}

    std::vector<double> IMoistureBC::R_Vector() const
    {
        const auto satOutside = saturationConcentrationAtTemperature(m_AirTemperature);
        const auto gconv = betaConv() * satOutside * m_AirHumidity;
        return m_PsiVector * gconv;
    }

    FenestrationCommon::SquareMatrix IMoistureBC::H_Matrix() const
    {
        std::vector<double> concentration(numOfBCNodes, 0);
        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            double T = m_Nodes[j].property(Variable::temperature);
            concentration[j] = saturationConcentrationAtTemperature(T);
        }
        const auto vaporFlux = concentration * betaConv();

        return m_PsiPsiMatrix.mmultRows(vaporFlux);
    }

    std::vector<double> IMoistureBC::betaConv() const
    {
        std::vector<double> betaCon(numOfBCNodes, 0);

        for(std::size_t j = 0; j < numOfBCNodes; ++j)
        {
            const double humidity = m_Nodes[j].property(Variable::humidity);
            if(humidity <= 1)
            {
                betaCon[j] = 1/(Constants::Cp_Air * Constants::Density_Air);
            }
            else
            {
                betaCon[j] = 0;
            }
        }
        return betaCon * convectiveCoefficient();
    }

    /////////////////////////////////////////////////////
    /// MoistureBCVariableHc
    /////////////////////////////////////////////////////
    MoistureBCVariableHc::MoistureBCVariableHc(size_t index1,
                                                   size_t index2,
                                                   const std::string & materialName,
                                                   double t_AirHumidity,
                                                   double t_AirTemperature) :
        IMoistureBC(index1, index2, materialName, t_AirHumidity, t_AirTemperature, ConvectionModel::Variable)
    {}

    std::vector<double> MoistureBCVariableHc::convectiveCoefficient() const
    {
        return m_ConvectiveCoeffCalc->value(m_Nodes, m_AirTemperature);
    }

    /////////////////////////////////////////////////////
    /// MoistureBCFixedHc
    /////////////////////////////////////////////////////
    MoistureBCFixedHc::MoistureBCFixedHc(size_t index1,
                                         size_t index2,
                                         const std::string & materialName,
                                         double t_AirHumidity,
                                         double t_AirTemperature,
                                         double convectiveCoefficient) :
        IMoistureBC(index1, index2, materialName, t_AirHumidity, t_AirTemperature, ConvectionModel::Fixed),
        m_ConvectiveCoefficient(convectiveCoefficient)
    {}

    std::vector<double> MoistureBCFixedHc::convectiveCoefficient() const
    {
        return m_ConvectiveCoeffCalc->value(m_Nodes, m_ConvectiveCoefficient);
    }
}   // namespace HygroThermFEM
