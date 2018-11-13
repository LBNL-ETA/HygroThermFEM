#include <cassert>
#include <iostream>

#include "Common.hxx"
#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "FEMunique.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM
{
    //////////////////////////////////////////////////////////////////////////////
    ///  IQLEMatrix2D
    //////////////////////////////////////////////////////////////////////////////
    IQLEIntegrator2D::IQLEIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        m_Global2D{t_Element},
        m_DPsiDxDyMatrix{numOfQuadrilateralNodes, SquareMatrix{numOfQuadrilateralNodes}}
    {}

    SquareMatrix IQLEIntegrator2D::integrate(const std::vector<double> & t_Values) const
    {
        const auto count = IntegrationPoints2D::Instance().count2D();

        SquareMatrix aMatrix(numOfQuadrilateralNodes);

        for(auto integrationPoint = 0u; integrationPoint < count; ++integrationPoint)
        {
            // Passing reference to matrix is faster then returning matrix from function
            // just because new matrix needs to be created and assigned. This was causing
            // some 40% slowdown.
            calculateMatrixInIntegrationPoint(aMatrix, t_Values, integrationPoint);
        }

        return aMatrix;
    }

    void IQLEIntegrator2D::calculateMatrixInIntegrationPoint(
      SquareMatrix & matrix,
      const std::vector<double> & t_Values,
      const std::size_t t_IntegrationPointIndex) const
    {
        assert(t_Values.size() == 4);
        assert(matrix.size() == 4);

        auto & intPointMatrix = m_DPsiDxDyMatrix[t_IntegrationPointIndex];

        for(size_t i = 0; i < t_Values.size(); ++i)
        {
            for(size_t j = 0; j < t_Values.size(); ++j)
            {
                matrix(i, j) += intPointMatrix(i, j) * 0.5 * (t_Values[i] + t_Values[j]);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDDUConductance2D
    //////////////////////////////////////////////////////////////////////////////

    QLEDDuIntegrator2D::QLEDDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            auto DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            auto DPsiDy = m_Global2D.DPsiDy(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto & DPsiDxDyMatrix = m_DPsiDxDyMatrix[integrationPoint];
            for(auto i = 0u; i < DPsiDxDyMatrix.size(); ++i)
            {
                for(auto j = 0u; j < DPsiDxDyMatrix.size(); ++j)
                {
                    DPsiDxDyMatrix(i, j) = (DPsiDx[i] * DPsiDx[j] + DPsiDy[i] * DPsiDy[j]) * det;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEConductanceDerivative2D
    //////////////////////////////////////////////////////////////////////////////

    QLEDpDuIntegrator2D::QLEDpDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {}

    void QLEDpDuIntegrator2D::setIndependentVariables(const std::vector<double> & t_Values)
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        auto & aElement = QuadrilateralLinearLocal2D::Instance();

        assert(t_Values.size() == numOfIntegrationPoints);

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto & psi = aElement.VPsi(integrationPoint);
            auto DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            auto DPsiDy = m_Global2D.DPsiDy(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto gammaX = 0.0;
            auto gammaY = 0.0;
            for(auto k = 0u; k < numOfIntegrationPoints; ++k)
            {
                gammaX += DPsiDx[k] * t_Values[k];
                gammaY += DPsiDy[k] * t_Values[k];
            }

            auto & psiPsiMatrix = m_DPsiDxDyMatrix[integrationPoint];
            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                for(auto j = 0u; j < numOfIntegrationPoints; ++j)
                {
                    psiPsiMatrix(i, j) =
                      det * (DPsiDx[i] * psi[j] * gammaX + DPsiDy[i] * psi[j] * gammaY);
                }
            }
        }
    }

    void QLEDpDuIntegrator2D::clearIntegrationMatrix()
    {
        m_DPsiDxDyMatrix.clear();
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLECapacitance2D
    //////////////////////////////////////////////////////////////////////////////

    QLECapacitanceIntegrator2D::QLECapacitanceIntegrator2D(
      const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        auto & aElement = QuadrilateralLinearLocal2D::Instance();

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto & psi = aElement.VPsi(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto & psiPsiMatrix = m_DPsiDxDyMatrix[integrationPoint];
            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                for(auto j = 0u; j < numOfIntegrationPoints; ++j)
                {
                    psiPsiMatrix(i, j) = det * psi[i] * psi[j];
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  DerivativeFunction
    //////////////////////////////////////////////////////////////////////////////

    DerivativeFunction::DerivativeFunction(iValue fixedValue, iValue derivativeValue) :
        fixedValue(std::move(fixedValue)),
        derivativeValue(std::move(derivativeValue))
    {}

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::IElementLinear2D(const size_t index1,
                                       const size_t index2,
                                       const size_t index3,
                                       const size_t index4,
                                       const std::string & materialName,
                                       const bool isLinear) :
        m_Material{MaterialPool::Instance().material(materialName)},
        m_Nodes{NodePool::Instance().getNode(index1),
                NodePool::Instance().getNode(index2),
                NodePool::Instance().getNode(index3),
                NodePool::Instance().getNode(index4)},
        m_Global2D{NodePool::Instance().getNode(index1),
                   NodePool::Instance().getNode(index2),
                   NodePool::Instance().getNode(index3),
                   NodePool::Instance().getNode(index4)},
        m_QLECapacitance2D{m_Global2D},
        m_DDuIntegrator{m_Global2D},
        m_Linear{isLinear}
    {
        auto matName = m_Material.name();

        /// Evaluating material influence in every node (This is important to know when
        /// calculating water content).
        while(!m_Nodes.last())
        {
            /// Form triangle of nodes. node1 is in center and angle is calculated at that node.
            auto & node1 = NodePool::Instance().getNode(m_Nodes.current().getNodeNumber());
            auto & node2 = NodePool::Instance().getNode(m_Nodes.previous().getNodeNumber());
            auto & node3 = NodePool::Instance().getNode(m_Nodes.next().getNodeNumber());

            /// Weighting coefficient depends on angle that is form by nodes next to node1.
            /// That coefficient is fraction of full circle.
            const auto weightingCoefficient =
              angleBetweenNodes(node1, node2, node3) / Constants::PI;
            /// Node will have possibility to calculate certain properties that will be
            /// material dependent.
            node1.assignMaterial(matName, weightingCoefficient);
            m_Nodes.moveToNext();
        }
    }

    FenestrationCommon::SquareMatrix IElementLinear2D::DDuMatrices() const
    {
        FenestrationCommon::SquareMatrix result{numOfQuadrilateralNodes};
        for(const auto & cond : m_DDuFunctions)
        {
            const auto values = cond->values(m_Nodes);
            result += m_DDuIntegrator.integrate(values);
        }

        return result;
    }

    FenestrationCommon::SquareMatrix IElementLinear2D::DpDuMatrices()
    {
        FenestrationCommon::SquareMatrix result{numOfQuadrilateralNodes};

        /// Integration matrix must be created every time because independent
        /// variables changed as well.

        auto count = 0u;
        m_QLEDpDuIntegrator2D.clear();
        for(const auto & cond : m_DpDuFunctions)
        {
            m_QLEDpDuIntegrator2D.emplace_back(m_Global2D);
            const auto aDerivatives = cond.derivativeValue->values(m_Nodes);
            m_QLEDpDuIntegrator2D[count].setIndependentVariables(aDerivatives);
            ++count;
        }

        /// Now rest of integration is performed as usual
        count = 0u;
        for(const auto & cond : m_DpDuFunctions)
        {
            const auto values = cond.fixedValue->values(m_Nodes);
            result += m_QLEDpDuIntegrator2D[count].integrate(values);
            ++count;
        }

        return result;
    }

    FenestrationCommon::SquareMatrix IElementLinear2D::capacitanceMatrices() const
    {
        FenestrationCommon::SquareMatrix result{numOfQuadrilateralNodes};
        for(const auto & cap : m_CapacitanceFunctions)
        {
            const auto values = cap->values(m_Nodes);
            result += m_QLECapacitance2D.integrate(values);
        }
        return result;
    }

    Node2D & IElementLinear2D::getNode(const std::size_t index)
    {
        assert(index < m_Nodes.size());
        return m_Nodes[index];
    }

    std::vector<std::size_t> IElementLinear2D::nodeIndexes() const
    {
        return m_Global2D.nodeIndexes();
    }

    const Material & IElementLinear2D::getMaterial() const
    {
        return m_Material;
    }

    bool IElementLinear2D::haveBothNodes( const size_t index1, const size_t index2 ) const
    {
        bool node1Found = false;
        bool node2Found = false;
        for(auto & node : m_Nodes)
        {
            node1Found = node1Found || node.get().getNodeNumber() == index1;
            node2Found = node2Found || node.get().getNodeNumber() == index2;
        }
        return node1Found && node2Found;
    }

    double IElementLinear2D::angleBetweenNodes(const Node2D & node1,
                                               const Node2D & node2,
                                               const Node2D & node3)
    {
        auto angle = std::abs(std::atan2(node3.Y() - node1.Y(), node3.X() - node1.X())
                              - std::atan2(node2.Y() - node1.Y(), node2.X() - node1.X()));
        if(angle > Constants::PI)
        {
            angle -= Constants::PI;
        }
        return angle;
    }

    bool IElementLinear2D::isLinear() const
    {
        return m_Linear;
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementThermalLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementThermalLinear2D::ElementThermalLinear2D( const size_t index1,
													const size_t index2,
													const size_t index3,
													const size_t index4,
													const std::string & materialName ) :
        IElementLinear2D(index1, index2, index3, index4, materialName)
    {
        //////////////////////////////////////////////////////////////////////////////////////
        /// Capacitance functions
        //////////////////////////////////////////////////////////////////////////////////////

        auto dryContent = (1 - m_Material.porosity()) * m_Material.density();
        StateValue liquidContent(Property::liquid);
        StateValue iceContent(Property::ice);
        // auto airContent = getMaterialAirFill(mat);
        StateValue airContent(Property::vapor);

        auto equivalentDensity =
          (dryContent * m_Material.density() + iceContent * Constants::Density_Ice
           + liquidContent * Constants::Density_Water + airContent * Constants::Density_Air)
          / (dryContent + iceContent + liquidContent + airContent);

        auto equivalentCapacitance =
          (dryContent * m_Material.heatCapacity() + iceContent * Constants::Cp_Ice
           + liquidContent * Constants::Cp_Water + airContent * Constants::Cp_Air)
          / (dryContent + iceContent + liquidContent + airContent);

        auto capacitance = equivalentDensity * equivalentCapacitance;

        m_CapacitanceFunctions.emplace_back(new decltype(capacitance)(capacitance));

        /// Conductance

        /// material
        auto materialConductivity = Constant(m_Material.thermalConductivity());

        /// vapor
        auto delta = Constant(2.5E-5 / m_Material.diffusionResistanceFactor());
        auto vaporConductivity = Constants::Cp_Vapor * delta * airContent;

        /// liquid
        auto humidity = StateValue(Property::humidity);
        auto liquidConductivity =
          SuctionCurve(m_Material.liquidTransportationCurve()) * Constants::Cp_Water * humidity;

        // iValue conductance = materialConductivity + vaporConductivity + liquidConductivity;
        auto conductance = materialConductivity + vaporConductivity + liquidConductivity;

        m_DDuFunctions.emplace_back(new decltype(conductance)(conductance));
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementMoistureLinear2D::ElementMoistureLinear2D( const size_t index1,
													  const size_t index2,
													  const size_t index3,
													  const size_t index4,
													  const std::string & materialName ) :
        IElementLinear2D(index1, index2, index3, index4, materialName, false)
    {
        //////////////////////////////////////////////////////////////////////////////
        /// Water vapor diffusion
        //////////////////////////////////////////////////////////////////////////////
        Constant delta(2.5E-5 / m_Material.diffusionResistanceFactor());
        auto conductance = delta * SaturationFunction();

        m_DDuFunctions.emplace_back(new decltype(conductance)(conductance));

        m_DpDuFunctions.emplace_back(std::unique_ptr<IValue>(new decltype(delta)(delta)),
                                     std::unique_ptr<IValue>(new SaturationFunction()));

        //////////////////////////////////////////////////////////////////////////////
        /// Water liquid transportation
        //////////////////////////////////////////////////////////////////////////////
        m_DDuFunctions.emplace_back(new SuctionCurve(m_Material.liquidTransportationCurve()));

        //////////////////////////////////////////////////////////////////////////////
        /// Creating capacitance function
        //////////////////////////////////////////////////////////////////////////////
        auto sorptionDerivative = TabularDerivative(m_Material.sorptionCurve(), Property::humidity);

        m_CapacitanceFunctions.emplace_back(new decltype(sorptionDerivative)(sorptionDerivative));
    }
}   // namespace MoisThermFEM
