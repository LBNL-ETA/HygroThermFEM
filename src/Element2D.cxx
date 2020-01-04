#include <cassert>

#include "Common.hxx"
#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "NodePool.hxx"
#include "MaterialPool.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"

namespace HygroThermFEM
{
    //////////////////////////////////////////////////////////////////////////////
    ///  IQLEMatrix2D
    //////////////////////////////////////////////////////////////////////////////
    IQLEIntegrator2D::IQLEIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        m_Global2D{t_Element},
        m_IntegrationMatrix{numOfQuadrilateralNodes, SquareMatrix{numOfQuadrilateralNodes}}
    {}

    SquareMatrix IQLEIntegrator2D::integrate(const std::vector<double> & t_Values) const
    {
        const auto count = IntegrationPoints2D::Instance().count2D();

        std::vector<std::vector<double>> aMatrix{numOfQuadrilateralNodes,
                                                 std::vector<double>(numOfQuadrilateralNodes, 0)};
        for(auto integrationPoint = 0u; integrationPoint < count; ++integrationPoint)
        {
            auto & intPointMatrix = m_IntegrationMatrix[integrationPoint];

            for(size_t i = 0; i < t_Values.size(); ++i)
            {
                for(size_t j = 0; j < t_Values.size(); ++j)
                {
                    aMatrix[i][j] += intPointMatrix(i, j) * 0.5 * (t_Values[i] + t_Values[j]);
                }
            }
        }

        return SquareMatrix{aMatrix};
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
            const auto & DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            const auto & DPsiDy = m_Global2D.DPsiDy(integrationPoint);
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
    ///  QLEDpDuIntegrator2D
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
            const auto & psi = aElement.Psi(integrationPoint);
            const auto & DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            const auto & DPsiDy = m_Global2D.DPsiDy(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto gammaX = 0.0;
            auto gammaY = 0.0;
            for(auto k = 0u; k < numOfIntegrationPoints; ++k)
            {
                gammaX += DPsiDx[k] * t_Values[k];
                gammaY += DPsiDy[k] * t_Values[k];
            }

            std::vector<std::vector<double>> matrix{numOfIntegrationPoints,
                                                    std::vector<double>(numOfIntegrationPoints, 0)};

            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                for(auto j = 0u; j < numOfIntegrationPoints; ++j)
                {
                    matrix[i][j] =
                      det * (psi[i] * DPsiDx[j] * gammaX + psi[i] * DPsiDy[j] * gammaY);
                }
            }
            m_IntegrationMatrix[integrationPoint] = SquareMatrix{matrix};
        }
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
            const auto & psi = aElement.Psi(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            for(auto i = 0u; i < numOfIntegrationPoints; ++i)
            {
                for(auto j = 0u; j < numOfIntegrationPoints; ++j)
                {
                    m_IntegrationMatrix[integrationPoint](i, j) = det * psi[i] * psi[j];
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::IElementLinear2D(const size_t index1,
                                       const size_t index2,
                                       const size_t index3,
                                       const size_t index4,
                                       const std::string & materialName,
                                       const Variable variable,
                                       const bool isLinear) :
        m_Material{MaterialPool::Instance().material(materialName)},
        m_FluxVariable(variable),
        m_Nodes{NodePool::Instance().getNode(index1),
                NodePool::Instance().getNode(index2),
                NodePool::Instance().getNode(index3),
                NodePool::Instance().getNode(index4)},
        m_Global2D{NodePool::Instance().getNode(index1),
                   NodePool::Instance().getNode(index2),
                   NodePool::Instance().getNode(index3),
                   NodePool::Instance().getNode(index4)},
        m_QLECapacitance2D{m_Global2D},
        m_Linear{isLinear && m_Material.isLinear()}
    {
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
              angleBetweenNodes(node1, node2, node3) / Constants::PI_HTF;
            /// Node will have possibility to calculate certain properties that will be
            /// material dependent.
            node1.assignMaterial(materialName, weightingCoefficient);
            m_Nodes.moveToNext();
        }
    }

    SquareMatrix IElementLinear2D::DDuMatrices() const
    {
        SquareMatrix result{numOfQuadrilateralNodes};

        const QLEDDuIntegrator2D DDuIntegrator{m_Global2D};
        for(const auto & cond : m_DDuFunctions)
        {
            const auto values = cond->values(m_Nodes);
            result += DDuIntegrator.integrate(values);
        }

        return result;
    }

    SquareMatrix IElementLinear2D::DpDuMatrices()
    {
        SquareMatrix result{numOfQuadrilateralNodes};

        /// Integration matrix must be created every time because independent
        /// variables changed as well.

        QLEDpDuIntegrator2D qleDpDuIntegrator2D{m_Global2D};
        for(const auto & cond : m_DpDuFunctions)
        {
            const auto aDerivatives = cond.derivativeValue->values(m_Nodes);
            qleDpDuIntegrator2D.setIndependentVariables(aDerivatives);
            const auto values = cond.fixedValue->values(m_Nodes);
            result += qleDpDuIntegrator2D.integrate(values);
        }

        return result;
    }

    SquareMatrix IElementLinear2D::capacitanceMatrices() const
    {
        SquareMatrix result{numOfQuadrilateralNodes};
        for(const auto & cap : m_CapacitanceFunctions)
        {
            const auto values = cap->values(m_Nodes);
            result += m_QLECapacitance2D.integrate(values);
        }
        return result;
    }

    std::vector<NodeFlux> IElementLinear2D::flux() const
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        std::vector<NodeFlux> results;
        for(const auto & cond : m_ConductanceFunctions)
        {
            const auto values = m_Nodes.properties(m_FluxVariable);
            const auto k = cond->values(m_Nodes);
            assert(k.size() == numOfQuadrilateralNodes);
            assert(values.size() == numOfQuadrilateralNodes);

            // First calculate Dt/Dx and Dt/Dy in Gauss points
            std::vector<double> VDtDx(numOfIntegrationPoints, 0);
            std::vector<double> VDtDy(numOfIntegrationPoints, 0);
            for(size_t i = 0; i < numOfIntegrationPoints; ++i)
            {
                const auto DPsiDx = m_Global2D.DPsiDx(i);
                const auto DPsiDy = m_Global2D.DPsiDy(i);
                assert(DPsiDx.size() == numOfIntegrationPoints);
                double DtDx{0};
                double DtDy{0};
                for(size_t j = 0; j < values.size(); ++j)
                {
                    DtDx += DPsiDx[j] * values[j];
                    DtDy += DPsiDy[j] * values[j];
                }
                VDtDx[i] = DtDx;
                VDtDy[i] = DtDy;
            }

            // Extrapolate flux from Gauss points to nodes.
            const auto a = 1.866025404;
            const auto b = -0.5;
            const auto c = 0.133974596;
            const std::vector<std::vector<double>> extrapolationCoefficients{
              {a, b, c, b}, {b, a, b, c}, {c, b, a, b}, {b, c, b, a}};

            for(size_t i = 0u; i < numOfQuadrilateralNodes; ++i)
            {
                auto valX{0.0};
                auto valY{0.0};
                for(const auto val : extrapolationCoefficients[i])
                {
                    valX -= k[i] * VDtDx[i] * val;
                    valY -= k[i] * VDtDy[i] * val;
                }
                results.emplace_back(valX, valY);
            }
        }
        return results;
    }

    Node2D & IElementLinear2D::getNode(const std::size_t index) const
    {
        assert(index < m_Nodes.size());
        return m_Nodes[index];
    }

    std::vector<std::size_t> IElementLinear2D::nodeIndexes() const
    {
        return m_Global2D.nodeIndexes();
    }

    const IMaterial & IElementLinear2D::getMaterial() const
    {
        return m_Material;
    }

    bool IElementLinear2D::haveBothNodes(const size_t index1, const size_t index2) const
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
                                               const Node2D & node3) const
    {
        auto angle = std::abs(std::atan2(node3.Y() - node1.Y(), node3.X() - node1.X())
                              - std::atan2(node2.Y() - node1.Y(), node2.X() - node1.X()));
        if(angle > Constants::PI_HTF)
        {
            angle -= Constants::PI_HTF;
        }
        return angle;
    }

    bool IElementLinear2D::isLinear() const
    {
        return m_Linear;
    }

    std::vector<double> IElementLinear2D::rightSideVector() const
    {
        std::vector<double> result(numOfQuadrilateralNodes, 0);

        const QLEDDuIntegrator2D DDuIntegrator{m_Global2D};

        /// SquareMatrix M{numOfQuadrilateralNodes};
        for(const auto & item : m_Matrix_x_Vector)
        {
            /// Calculate functions base on node properties
            const auto values = item.MatrixFunction->values(m_Nodes);
            /// And then integrate them
            auto M = DDuIntegrator.integrate(values);
            auto B = m_Nodes.properties(item.PropertyVector);

            result = result + M * B;
        }

        return result;
    }

    IElementLinear2D::MatrixVector::MatrixVector(iValue && matrixFunction,
                                                 const Variable propertyVector) :
        MatrixFunction(std::move(matrixFunction)), PropertyVector(propertyVector)
    {}

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D::DerivativeFunction
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::DerivativeFunction::DerivativeFunction(iValue fixedValue,
                                                             iValue derivativeValue) :
        fixedValue(std::move(fixedValue)), derivativeValue(std::move(derivativeValue))
    {}

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementThermalLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementThermalLinear2D::ElementThermalLinear2D(const size_t index1,
                                                   const size_t index2,
                                                   const size_t index3,
                                                   const size_t index4,
                                                   const std::string & materialName) :
        IElementLinear2D(index1, index2, index3, index4, materialName, Variable::temperature)
    {
        //////////////////////////////////////////////////////////////////////////////////////
        /// Capacitance functions
        //////////////////////////////////////////////////////////////////////////////////////

        // const auto dryContent = (1 - m_Material.porosity()) * m_Material.density();
        const StateValue liquidContent(Variable::liquid);
        const StateValue iceContent(Variable::ice);
        // auto airContent = getMaterialAirFill(mat);
        // const StateValue airContent(Variable::vapor);

        // Vapor content changes are ignored for now
        const auto equivalentDensity = m_Material.density() + liquidContent + iceContent;

        const auto equivalentCapacitance =
          m_Material.heatCapacity() + (iceContent / Constants::Density_Ice) * Constants::Cp_Ice
          + (liquidContent / Constants::Density_Water) * Constants::Cp_Water;

        auto capacitance = equivalentDensity * equivalentCapacitance;

        Cap(capacitance);

        // Phase change part
        // This is incorrect phase change equation. Correct one is kept in branch IceContentFix.
        // Disable this for now because solver is not producing correct results
        // auto waterWithoutVapor = liquidContent + iceContent;
        // Cap(PhaseChange() * waterWithoutVapor);

        //////////////////////////////////////////////////////////////////////////
        /// Conductance
        //////////////////////////////////////////////////////////////////////////

        /// material
        // const auto materialConductivity =
        //  TabularFunction1D(m_Material.thermalConductivityMoistureAndTemperatureDependent(), Variable::water);

        if(SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent())
        {
            auto materialConductivity = m_Material.thermalConductivityMoistureAndTemperatureDependent();
            DDu(materialConductivity);
        }
        else
        {
            auto materialConductivity = Constant(m_Material.thermalConductivityDry());
            DDu(materialConductivity);
        }


        //////////////////////////////////////////////////////////////////////
        ///  Conversion from liquid to gas (vapor part)
        //////////////////////////////////////////////////////////////////////
        const auto delta = Constant(2.5E-5 / m_Material.diffusionResistanceFactor());
        if(!SimulationProperties::Instance().excludeHeatOfEvaporation())
        {
            auto h = HeatOfEvaporation() * delta;

            multiplies(h, Variable::vapor);
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conversion from liquid to gas (air part)
        //////////////////////////////////////////////////////////////////////

        /// TODO: Add this later when air pressure equation is added
        // auto waterVaporPressure = SaturationFunction() * StateValue(Variable::humidity);

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from liquid
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeCapillaryConduction())
        {
            auto humidity = StateValue(Variable::humidity);
            const TabularDerivativeSmooth sorptionDerivative(m_Material.sorptionCurve(),
                                                             Variable::humidity);
            const LiquidTransportationCurve Dl(m_Material.liquidTransportationCurve());
            auto cd = Constant(-1) * Dl * sorptionDerivative * Constants::Cp_Water;
            DpDu(cd, humidity);
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from vapor
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeVaporDiffusionConduction())
        {
            auto vapCond = Constant(-1) * delta * Constants::Cp_Vapor;
            StateValue vaporContent(Variable::vapor);

            DpDu(vapCond, vaporContent);
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from airflow
        //////////////////////////////////////////////////////////////////////

        //////////////////////////////////////////////////////////////////////
        /// SolidMaterial conductance for flux calculations
        //////////////////////////////////////////////////////////////////////
        if(SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent())
        {
            Constant matCond{m_Material.thermalConductivityDry()};
            CondFlux(matCond);
        }
        else
        {
            auto matCond{m_Material.thermalConductivityMoistureAndTemperatureDependent()};
            CondFlux(matCond);
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementMoistureLinear2D::ElementMoistureLinear2D(const size_t index1,
                                                     const size_t index2,
                                                     const size_t index3,
                                                     const size_t index4,
                                                     const std::string & materialName) :
        IElementLinear2D(index1, index2, index3, index4, materialName, Variable::humidity, false)
    {
        //////////////////////////////////////////////////////////////////////////////
        /// Water vapor diffusion
        //////////////////////////////////////////////////////////////////////////////
        Constant delta(2.5E-5 / m_Material.diffusionResistanceFactor());
        auto conductance = delta * SaturationFunction();

        // Note that diffusion equation is partial derivative that gets split into two derivative
        // terms which then get into system as DDU and DpDu parts.
        DDu(conductance);

        DpDu(delta, SaturationFunction());

        //////////////////////////////////////////////////////////////////////////////
        /// Water liquid transportation
        //////////////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeWaterLiquidTransportation())
        {
            auto sorptionDerivative1 =
              TabularDerivativeSmooth(m_Material.sorptionCurve(), Variable::humidity);
            auto WaterLiquidTransport =
              LiquidTransportationCurve(m_Material.liquidTransportationCurve())
              * sorptionDerivative1;
            DDu(WaterLiquidTransport);
        }

        //////////////////////////////////////////////////////////////////////////////
        /// Creating capacitance function
        //////////////////////////////////////////////////////////////////////////////
        auto sorptionDerivative2 =
          TabularDerivativeSmooth(m_Material.sorptionCurve(), Variable::humidity);
        Cap(sorptionDerivative2);

        //////////////////////////////////////////////////////////////////////
        /// Functions for flux calculations
        //////////////////////////////////////////////////////////////////////
        CondFlux(conductance);
        // Cond(delta);
        CondFlux(LiquidTransportationCurve(m_Material.liquidTransportationCurve()));
    }

}   // namespace HygroThermFEM
