#include <cassert>
#include <iostream>

#include "Common.hxx"
#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "MaterialProperties.hxx"
#include "NodePool.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "FEMunique.hxx"

using FenestrationCommon::SquareMatrix;

namespace MoisThermFEM
{
    //////////////////////////////////////////////////////////////////////////////
    ///  IQLEMatrix2D
    //////////////////////////////////////////////////////////////////////////////
    IQLEMatrix2D::IQLEMatrix2D(const QuadrilateralLinearGlobal2D & t_Element) :
        m_Global2D{t_Element},
        m_IntegrationMatrix{numOfQuadrilateralNodes, SquareMatrix{numOfQuadrilateralNodes}}
    {}

    SquareMatrix IQLEMatrix2D::integrate(const std::vector<double> & t_Values) const
    {
        const auto count = IntegrationPoints2D::Instance().count2D();

        // SquareMatrix aMatrix{ numOfQuadrilateralNodes };
        std::vector<std::vector<double>> aMatrix(numOfQuadrilateralNodes,
                                                 std::vector<double>(numOfQuadrilateralNodes, 0));

        for(auto i = 0u; i < count; ++i)
        {
            calculateMatrixInIntegrationPoint(t_Values, i, aMatrix);
        }

        return SquareMatrix{aMatrix};
    }

    void IQLEMatrix2D::calculateMatrixInIntegrationPoint(
      const std::vector<double> & t_Values,
      const std::size_t t_IntegrationPointIndex,
      std::vector<std::vector<double>> & t_Matrix) const
    {
        assert(t_Values.size() == 4);

        auto & intPointMatrix = m_IntegrationMatrix[t_IntegrationPointIndex];

        for(size_t i = 0; i < t_Matrix.size(); ++i)
        {
            for(size_t j = 0; j < t_Matrix.size(); ++j)
            {
                t_Matrix[i][j] += intPointMatrix(i, j) * 0.5 * (t_Values[i] + t_Values[j]);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEConductance2D
    //////////////////////////////////////////////////////////////////////////////

    QLEConductance2D::QLEConductance2D(const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEMatrix2D{t_Element}
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            auto DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            auto DPsiDy = m_Global2D.DPsiDy(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto & DPsiDxDyMatrix = m_IntegrationMatrix[integrationPoint];
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

    QLEConductanceDerivative2D::QLEConductanceDerivative2D(
      const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEMatrix2D{t_Element}
    {}

    void QLEConductanceDerivative2D::updateIntegrationMatrix(const std::vector<double> & t_Values)
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

            auto & psiPsiMatrix = m_IntegrationMatrix[integrationPoint];
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

    void QLEConductanceDerivative2D::clearIntegrationMatrix()
    {
        m_IntegrationMatrix.clear();
    }

	//////////////////////////////////////////////////////////////////////////////
    ///  QLECapacitance2D
    //////////////////////////////////////////////////////////////////////////////

    QLECapacitance2D::QLECapacitance2D(const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEMatrix2D{t_Element}
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        auto & aElement = QuadrilateralLinearLocal2D::Instance();

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto & psi = aElement.VPsi(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto & psiPsiMatrix = m_IntegrationMatrix[integrationPoint];
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

    DerivativeFunction::DerivativeFunction(iValue fixedTerm, iValue derivativeTerm) :
        fixedTerm(std::move(fixedTerm)),
        derivativeTerm(std::move(derivativeTerm))
    {}

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::IElementLinear2D(const Node2D & t_Node1,
                                       const Node2D & t_Node2,
                                       const Node2D & t_Node3,
                                       const Node2D & t_Node4,
                                       const Material & t_Material) :
        m_Material{t_Material},
        m_Node{{t_Node1, t_Node2, t_Node3, t_Node4}},
        m_Global2D{t_Node1, t_Node2, t_Node3, t_Node4},
        m_QLECapacitance2D{m_Global2D},
        m_QLEConductance2D{m_Global2D}
    {
        auto matName = m_Material.name();

        /// Iterating through unique nodes in element. Note that element can be triangular in
        /// which case one of the node numbers will be repeated twice.
        while(!m_Node.last())
        {
            auto & node1 = NodePool::Instance().getNode(m_Node.current().getNodeNumber());
            auto & node2 = NodePool::Instance().getNode(m_Node.previous().getNodeNumber());
            auto & node3 = NodePool::Instance().getNode(m_Node.next().getNodeNumber());

            /// Weighting coefficient depends on angle that is form by nodes next to node1.
            /// That coefficient is fraction of full circle.
            const auto weightingCoefficient =
              angleBetweenNodes(node1, node2, node3) / Constants::PI;
            /// Node will have possibility to calculate certain properties that will be
            /// material dependent.
            node1.assignMaterial(matName, weightingCoefficient);
            m_Node.moveToNext();
        }
    }

    FenestrationCommon::SquareMatrix IElementLinear2D::conductanceMatrix() const
    {
        FenestrationCommon::SquareMatrix result{numOfQuadrilateralNodes};
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        for(const auto & cond : m_Conductance)
        {
            std::vector<double> values(numOfIntegrationPoints);
            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                values[i] = cond->value(m_Node[i].getState());
            }
            result += m_QLEConductance2D.integrate(values);
        }

        return result;
    }

    FenestrationCommon::SquareMatrix IElementLinear2D::conductanceDerivativeMatrix()
    {
        FenestrationCommon::SquareMatrix result{numOfQuadrilateralNodes};
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();

        /// Integration matrix must be created every time
        std::vector<double> aDerivatives(numOfIntegrationPoints);

        auto count = 0u;
        m_QLEDerivativeConductance.clear();
        for(const auto & cond : m_DerivativeConductance)
        {
            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                aDerivatives[i] = cond.derivativeTerm->value(m_Node[i].getState());
            }
            m_QLEDerivativeConductance.emplace_back(m_Global2D);
            m_QLEDerivativeConductance[count].updateIntegrationMatrix(aDerivatives);
            ++count;
        }

        /// Now rest of integration is performed as usual
        count = 0u;
        for(const auto & cond : m_DerivativeConductance)
        {
            std::vector<double> values(numOfIntegrationPoints);
            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                values[i] = cond.fixedTerm->value(m_Node[i].getState());
            }
            result += m_QLEDerivativeConductance[count].integrate(values);
            ++count;
        }

        return result;
    }

    FenestrationCommon::SquareMatrix IElementLinear2D::capacitanceMatrix() const
    {
        FenestrationCommon::SquareMatrix result{numOfQuadrilateralNodes};
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        for(const auto & cap : m_Capacitance)
        {
            std::vector<double> values(numOfIntegrationPoints);
            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                values[i] = cap->value(m_Node[i].getState());
            }
            result += m_QLECapacitance2D.integrate(values);
        }
        return result;
    }

    Node2D & IElementLinear2D::getNode(const std::size_t index)
    {
        assert(index < m_Node.size());
        return m_Node[index];
    }

    std::vector<std::size_t> IElementLinear2D::nodeIndexes() const
    {
        return m_Global2D.nodeIndexes();
    }

    const Material & IElementLinear2D::getMaterial() const
    {
        return m_Material;
    }

    bool IElementLinear2D::haveBothNodes(const Node2D & t_Node1, const Node2D & t_Node2) const
    {
        bool node1Found = false;
        bool node2Found = false;
        for(auto & node : m_Node)
        {
            node1Found = node1Found || node == t_Node1;
            node2Found = node2Found || node == t_Node2;
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

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementThermalLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementThermalLinear2D::ElementThermalLinear2D(const Node2D & t_Node1,
                                                   const Node2D & t_Node2,
                                                   const Node2D & t_Node3,
                                                   const Node2D & t_Node4,
                                                   const Material & mat) :
        IElementLinear2D(t_Node1, t_Node2, t_Node3, t_Node4, mat)
    {
        auto waterFill = getLiquidWaterFill(mat);
        auto airFill = getAirFill(mat);

        auto waterCapacitance = waterFill * (Constants::Density_Water * Constants::Cp_Water);
        auto airCapacitance = airFill * (Constants::Density_Air * Constants::Cp_Air);
        auto dryCapacitance = (1 - mat.porosity()) * (mat.density() * mat.heatCapacity());

        auto capacitance = waterCapacitance + airCapacitance + dryCapacitance;

        m_Capacitance.emplace_back(new decltype(capacitance)(capacitance));

        /// Conductance

        /// material
        auto materialConductivity = Constant(mat.thermalConductivity());

        /// vapor
        auto delta = Constant(2.5E-5 / mat.diffusionResistanceFactor());
        auto vaporConductivity = Constants::Cp_Vapor * delta * airFill;

        /// liquid
        auto humidity = StateValue(Property::humidity);
        auto liquidConductivity =
          SuctionFunction(mat.liquidTransportationCurve(), Property::humidity) * Constants::Cp_Water
          * humidity;

        // iValue conductance = materialConductivity + vaporConductivity + liquidConductivity;
        auto conductance = materialConductivity + vaporConductivity + liquidConductivity;

        m_Conductance.emplace_back(new decltype(conductance)(conductance));
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementMoistureLinear2D::ElementMoistureLinear2D(const Node2D & t_Node1,
                                                     const Node2D & t_Node2,
                                                     const Node2D & t_Node3,
                                                     const Node2D & t_Node4,
                                                     const Material & mat) :
        IElementLinear2D(t_Node1, t_Node2, t_Node3, t_Node4, mat)
    {
        //////////////////////////////////////////////////////////////////////////////
        /// Creating conductance function for vapor
        //////////////////////////////////////////////////////////////////////////////
        auto delta = Constant(2.5E-5 / mat.diffusionResistanceFactor());
        auto saturationFunction = SaturationFunction(Property::temperature);
        auto conductance = delta * saturationFunction;

        m_Conductance.emplace_back(new decltype(conductance)(conductance));

        m_DerivativeConductance.emplace_back(
          std::unique_ptr<IValue>(new decltype(delta)(delta)),
          std::unique_ptr<IValue>(new decltype(saturationFunction)(saturationFunction)));

        //////////////////////////////////////////////////////////////////////////////
        /// Creating conductance function for liquid
        //////////////////////////////////////////////////////////////////////////////
        auto suctionCurve = SuctionFunction(mat.liquidTransportationCurve(), Property::humidity);
        m_Conductance.emplace_back(new decltype(suctionCurve)(suctionCurve));

        //////////////////////////////////////////////////////////////////////////////
        /// Creating capacitance function
        //////////////////////////////////////////////////////////////////////////////
        auto sorptionDerivative = TabularDerivative(mat.sorptionCurve(), Property::humidity);

        m_Capacitance.emplace_back(new decltype(sorptionDerivative)(sorptionDerivative));
    }
}   // namespace MoisThermFEM
