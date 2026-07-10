#include <cassert>
#include <numbers>

#include "Common.hxx"
#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "Nodes.hxx"
#include "Materials.hxx"
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

            for(size_t row = 0; row < t_Values.size(); ++row)
            {
                for(size_t col = 0; col < t_Values.size(); ++col)
                {
                    aMatrix[row][col] += intPointMatrix(row, col)
                                         * 0.5 * (t_Values[row] + t_Values[col]);
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
            for(auto row = 0u; row < DPsiDxDyMatrix.size(); ++row)
            {
                for(auto col = 0u; col < DPsiDxDyMatrix.size(); ++col)
                {
                    DPsiDxDyMatrix(row, col) =
                      (DPsiDx[row] * DPsiDx[col] + DPsiDy[row] * DPsiDy[col]) * det;
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
            for(auto idx = 0u; idx < numOfIntegrationPoints; ++idx)
            {
                gammaX += DPsiDx[idx] * t_Values[idx];
                gammaY += DPsiDy[idx] * t_Values[idx];
            }

            std::vector<std::vector<double>> matrix{numOfIntegrationPoints,
                                                    std::vector<double>(numOfIntegrationPoints, 0)};

            for(auto row = 0u; row < numOfIntegrationPoints; ++row)
            {
                for(auto col = 0u; col < numOfIntegrationPoints; ++col)
                {
                    matrix[row][col] =
                      det * (psi[row] * DPsiDx[col] * gammaX + psi[row] * DPsiDy[col] * gammaY);
                }
            }
            m_IntegrationMatrix[integrationPoint] = SquareMatrix{matrix};
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDpDuConsistentIntegrator2D
    //////////////////////////////////////////////////////////////////////////////

    QLEDpDuConsistentIntegrator2D::QLEDpDuConsistentIntegrator2D(
      const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {}

    void QLEDpDuConsistentIntegrator2D::setIndependentVariables(
      const std::vector<double> & t_Values)
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
            for(auto idx = 0u; idx < numOfIntegrationPoints; ++idx)
            {
                gammaX += DPsiDx[idx] * t_Values[idx];
                gammaY += DPsiDy[idx] * t_Values[idx];
            }

            std::vector<std::vector<double>> matrix{numOfIntegrationPoints,
                                                    std::vector<double>(numOfIntegrationPoints, 0)};

            for(auto row = 0u; row < numOfIntegrationPoints; ++row)
            {
                for(auto col = 0u; col < numOfIntegrationPoints; ++col)
                {
                    // Transpose of QLEDpDuIntegrator2D: the test function (row) is the one
                    // that is differentiated -> integral( (grad psi_row . grad p) psi_col ).
                    matrix[row][col] =
                      det * (DPsiDx[row] * gammaX + DPsiDy[row] * gammaY) * psi[col];
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

            for(auto row = 0u; row < numOfIntegrationPoints; ++row)
            {
                for(auto col = 0u; col < numOfIntegrationPoints; ++col)
                {
                    m_IntegrationMatrix[integrationPoint](row, col) =
                      det * psi[row] * psi[col];
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::IElementLinear2D(Nodes & nodePool,
                                       Materials & materialPool,
                                       const size_t index1,
                                       const size_t index2,
                                       const size_t index3,
                                       const size_t index4,
                                       const std::string & materialName,
                                       const Variable variable,
                                       const bool isLinear) :
        m_Material{materialPool.material(materialName)},
        m_FluxVariable(variable),
        m_Nodes{nodePool.getNode(index1),
                nodePool.getNode(index2),
                nodePool.getNode(index3),
                nodePool.getNode(index4)},
        m_Global2D{nodePool.getNode(index1),
                   nodePool.getNode(index2),
                   nodePool.getNode(index3),
                   nodePool.getNode(index4)},
        m_QLECapacitance2D{m_Global2D},
        m_Linear{isLinear && m_Material.isLinear()}
    {
        /// Evaluating material influence in every node (This is important to know when
        /// calculating water content).
        // TODO: This will perform for every domain. Consequence is that weighting factors will be
        // doubled when two subdomains are used or tripled if three domains are used. However,
        // results will be correct since ratio will be correct.
        while(!m_Nodes.last())
        {
            /// Form triangle of nodes. node1 is in center and angle is calculated at that node.
            auto & node1 = nodePool.getNode(m_Nodes.current().getNodeNumber());
            auto & node2 = nodePool.getNode(m_Nodes.previous().getNodeNumber());
            auto & node3 = nodePool.getNode(m_Nodes.next().getNodeNumber());

            /// Weighting coefficient depends on angle that is form by nodes next to node1.
            /// That coefficient is fraction of full circle.
            const auto weightingCoefficient =
              angleBetweenNodes(node1, node2, node3) / (2 * std::numbers::pi);
            /// Node will have possibility to calculate certain properties that will be
            /// material dependent. Pass the already-resolved material reference directly.
            node1.assignMaterial(m_Material, weightingCoefficient);
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

    SquareMatrix IElementLinear2D::DpDuMatrices() const
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

        // Consistent (integrated-by-parts) coupling terms. Only the moisture element
        // registers these (for the vapour temperature-gradient term); thermal elements
        // leave m_DpDuConsistentFunctions empty and are unaffected. See D1.
        QLEDpDuConsistentIntegrator2D qleDpDuConsistentIntegrator2D{m_Global2D};
        for(const auto & cond : m_DpDuConsistentFunctions)
        {
            const auto aDerivatives = cond.derivativeValue->values(m_Nodes);
            qleDpDuConsistentIntegrator2D.setIndependentVariables(aDerivatives);
            const auto values = cond.fixedValue->values(m_Nodes);
            result += qleDpDuConsistentIntegrator2D.integrate(values);
        }

        return result;
    }

    SquareMatrix IElementLinear2D::capacitanceMatrices() const
    {
        SquareMatrix result{numOfQuadrilateralNodes};
        for(const auto & cap : m_CapacitanceFunctions)
        {
            const auto values = cap->values(m_Nodes);
            result += m_LumpCapacityNodally ? nodalLumpedCapacity(values)
                                            : m_QLECapacitance2D.integrate(values);
        }
        return result;
    }

    SquareMatrix
      IElementLinear2D::nodalLumpedCapacity(const std::vector<double> & nodalCapacity) const
    {
        // Nodal volumes v_i = integral(psi_i) = row sums of the consistent mass matrix
        // integral(psi_i psi_j) (obtained by integrating with a unit coefficient).
        const auto massMatrix =
          m_QLECapacitance2D.integrate(std::vector<double>(numOfQuadrilateralNodes, 1.0));

        SquareMatrix result{numOfQuadrilateralNodes};
        for(std::size_t i = 0; i < numOfQuadrilateralNodes; ++i)
        {
            double nodalVolume = 0.0;
            for(std::size_t j = 0; j < numOfQuadrilateralNodes; ++j)
            {
                nodalVolume += massMatrix(i, j);
            }
            result(i, i) = nodalVolume * nodalCapacity[i];
        }
        return result;
    }

    std::vector<std::array<double, 2>> IElementLinear2D::stateGradientsAtGaussPoints() const
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        const auto values = m_Nodes.properties(m_FluxVariable);
        assert(values.size() == numOfQuadrilateralNodes);

        std::vector<std::array<double, 2>> gradients(numOfIntegrationPoints, {0.0, 0.0});
        for(size_t gpIdx = 0; gpIdx < numOfIntegrationPoints; ++gpIdx)
        {
            const auto DPsiDx = m_Global2D.DPsiDx(gpIdx);
            const auto DPsiDy = m_Global2D.DPsiDy(gpIdx);
            assert(DPsiDx.size() == numOfIntegrationPoints);
            double DtDx{0};
            double DtDy{0};
            for(size_t nodeIdx = 0; nodeIdx < values.size(); ++nodeIdx)
            {
                DtDx += DPsiDx[nodeIdx] * values[nodeIdx];
                DtDy += DPsiDy[nodeIdx] * values[nodeIdx];
            }
            gradients[gpIdx] = {DtDx, DtDy};
        }
        return gradients;
    }

    std::vector<NodeFlux> IElementLinear2D::flux() const
    {
        const auto gradients = stateGradientsAtGaussPoints();
        std::vector<NodeFlux> results;
        for(const auto & cond : m_ConductanceFunctions)
        {
            const auto conductivityValues = cond->values(m_Nodes);
            assert(conductivityValues.size() == numOfQuadrilateralNodes);

            // Extrapolate flux from 2x2 Gauss points to corner nodes using the
            // standard bilinear extrapolation matrix.  Coefficients derive from
            // evaluating the Gauss-point shape functions at the corner nodes:
            //   A = 1 + sqrt(3)/2,  B = -1/2,  C = 1 - sqrt(3)/2
            constexpr auto gaussToNodeA = 1.866025404;
            constexpr auto gaussToNodeB = -0.5;
            constexpr auto gaussToNodeC = 0.133974596;
            const std::vector<std::vector<double>> extrapolationCoefficients{
              {gaussToNodeA, gaussToNodeB, gaussToNodeC, gaussToNodeB},
              {gaussToNodeB, gaussToNodeA, gaussToNodeB, gaussToNodeC},
              {gaussToNodeC, gaussToNodeB, gaussToNodeA, gaussToNodeB},
              {gaussToNodeB, gaussToNodeC, gaussToNodeB, gaussToNodeA}};

            for(size_t nodeIdx = 0u; nodeIdx < numOfQuadrilateralNodes; ++nodeIdx)
            {
                auto fluxX{0.0};
                auto fluxY{0.0};
                for(const auto coeff : extrapolationCoefficients[nodeIdx])
                {
                    fluxX -= conductivityValues[nodeIdx] * gradients[nodeIdx][0] * coeff;
                    fluxY -= conductivityValues[nodeIdx] * gradients[nodeIdx][1] * coeff;
                }
                results.emplace_back(fluxX, fluxY);
            }
        }
        return results;
    }

    std::vector<NodeFlux> IElementLinear2D::fluxAtGaussPoints() const
    {
        const auto gradients = stateGradientsAtGaussPoints();
        const auto cond = meanConductivity();
        std::vector<NodeFlux> result;
        result.reserve(gradients.size());
        for(const auto & grad : gradients)
        {
            result.emplace_back(-cond * grad[0], -cond * grad[1]);
        }
        return result;
    }

    std::vector<std::array<double, 2>> IElementLinear2D::gaussPointGlobalCoordinates() const
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        std::vector<std::array<double, 2>> coords;
        coords.reserve(numOfIntegrationPoints);
        for(size_t gpIdx = 0; gpIdx < numOfIntegrationPoints; ++gpIdx)
        {
            coords.push_back({m_Global2D.xg(gpIdx), m_Global2D.yg(gpIdx)});
        }
        return coords;
    }

    double IElementLinear2D::meanConductivity() const
    {
        double total{0.0};
        std::size_t count{0};
        for(const auto & cond : m_ConductanceFunctions)
        {
            for(const auto value : cond->values(m_Nodes))
            {
                total += value;
                ++count;
            }
        }
        return count == 0 ? 0.0 : total / static_cast<double>(count);
    }

    INode2D & IElementLinear2D::getNode(const std::size_t index) const
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
        if(angle > std::numbers::pi)
        {
            angle -= std::numbers::pi;
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

    ElementThermalLinear2D::ElementThermalLinear2D(Nodes & nodePool,
                                                   Materials & materialPool,
                                                   const size_t index1,
                                                   const size_t index2,
                                                   const size_t index3,
                                                   const size_t index4,
                                                   const std::string & materialName) :
        IElementLinear2D(nodePool, materialPool, index1, index2, index3, index4, materialName, Variable::temperature)
    {
        //////////////////////////////////////////////////////////////////////////////////////
        /// Capacitance functions
        //////////////////////////////////////////////////////////////////////////////////////

        if(m_Material.hasDensity())
        {
            const StateValue liquidContent(Variable::liquid);
            const StateValue iceContent(Variable::ice);

            // Vapor content changes are ignored for now
            const auto equivalentDensity = m_Material.density() + liquidContent + iceContent;

            const auto equivalentCapacitance =
              m_Material.heatCapacity() + (iceContent / Constants::Density_Ice) * Constants::Cp_Ice
              + (liquidContent / Constants::Density_Water) * Constants::Cp_Water;

            auto capacitance = equivalentDensity * equivalentCapacitance;

            Cap(capacitance);
        }

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
        //  TabularFunction1D(m_Material.thermalConductivityMoistureAndTemperatureDependent(),
        //  Variable::water);

        if(SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent()
           && m_Material.hasThermalConductivityMoistureAndTemperatureDependent())
        {
            auto materialConductivity =
              m_Material.thermalConductivityMoistureAndTemperatureDependent();
            DDu(materialConductivity);

            auto matCond{m_Material.thermalConductivityMoistureAndTemperatureDependent()};
            CondFlux(matCond);
        }
        if(!SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent()
           && m_Material.hasThermalConductivityDry())
        {
            auto materialConductivity = Constant(m_Material.thermalConductivityDry());
            DDu(materialConductivity);

            Constant matCond{m_Material.thermalConductivityDry()};
            CondFlux(matCond);
        }


        //////////////////////////////////////////////////////////////////////
        ///  Conversion from liquid to gas (vapor part)
        //////////////////////////////////////////////////////////////////////
        if(m_Material.hasDiffusionResistanceFactor())
        {
            const auto delta = Constant(2.5E-5 / m_Material.diffusionResistanceFactor());
            if(!SimulationProperties::Instance().excludeHeatOfEvaporation())
            {
                auto evaporationHeat = HeatOfEvaporation() * delta;

                multiplies(evaporationHeat, Variable::vapor);
            }
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conversion from liquid to gas (air part)
        //////////////////////////////////////////////////////////////////////

        /// TODO: Add this later when air pressure equation is added
        // auto waterVaporPressure = SaturationFunction() * StateValue(Variable::humidity);

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from liquid
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeCapillaryConduction()
           && m_Material.hasSorptionCurve() && m_Material.hasLiquidTransportationCurve())
        {
            auto humidity = StateValue(Variable::humidity);
            const TabularDerivativeSmooth sorptionDerivative(m_Material.sorptionCurve(),
                                                             Variable::humidity);
            const LiquidTransportationCurve Dl(m_Material.liquidTransportationCurve(), m_Material);
            auto cd = Constant(-1) * Dl * sorptionDerivative * Constants::Cp_Water;
            DpDu(cd, humidity);
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from vapor
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeVaporDiffusionConduction()
           && m_Material.hasDiffusionResistanceFactor())
        {
            const auto delta = Constant(2.5E-5 / m_Material.diffusionResistanceFactor());
            auto vapCond = Constant(-1) * delta * Constants::Cp_Vapor;
            StateValue vaporContent(Variable::vapor);

            DpDu(vapCond, vaporContent);
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementMoistureLinear2D::ElementMoistureLinear2D(Nodes & nodePool,
                                                     Materials & materialPool,
                                                     const size_t index1,
                                                     const size_t index2,
                                                     const size_t index3,
                                                     const size_t index4,
                                                     const std::string & materialName) :
        IElementLinear2D(nodePool, materialPool, index1, index2, index3, index4, materialName, Variable::humidity, false)
    {
        // Lump the moisture capacity nodally (diag(v_i * xi_i)) so the secant capacity
        // conserves stored moisture exactly; the averaged-then-lumped form only conserves
        // approximately and leaks ~1% under strong temperature gradients. See D6-refinement.
        m_LumpCapacityNodally = true;

        //////////////////////////////////////////////////////////////////////////////
        /// Water vapor diffusion
        //////////////////////////////////////////////////////////////////////////////
        if(m_Material.hasDiffusionResistanceFactor())
        {
            Constant delta(2.5E-5 / m_Material.diffusionResistanceFactor());
            auto conductance = delta * SaturationFunction();

            // The vapour term grad.(delta grad(phi c_sat)) splits into a moisture-gradient
            // half (delta c_sat grad phi) and a temperature-gradient half (delta phi grad
            // c_sat). Both must be integrated by parts consistently: the first is the DDu
            // (stiffness) term, the second is the *consistent* DpDu term (test function
            // differentiated). Using the plain DpDu here assembled the transpose and dropped
            // a term, which made results swing with the temperature field. See D1.
            DDu(conductance);

            DpDuConsistent(delta, SaturationFunction());

            // Function for flux calculations.
            CondFlux(conductance);
        }

        //////////////////////////////////////////////////////////////////////////////
        /// Water liquid transportation
        //////////////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeWaterLiquidTransportation()
           && m_Material.hasSorptionCurve() && m_Material.hasLiquidTransportationCurve())
        {
            auto sorptionDerivative =
              TabularDerivativeSmooth(m_Material.sorptionCurve(), Variable::humidity);
            auto WaterLiquidTransport =
              LiquidTransportationCurve(m_Material.liquidTransportationCurve(), m_Material)
              * sorptionDerivative;
            DDu(WaterLiquidTransport);
        }

        //////////////////////////////////////////////////////////////////////////////
        /// Creating capacitance function
        //////////////////////////////////////////////////////////////////////////////
        if(m_Material.hasSorptionCurve())
        {
            // Mass-conservative (secant) storage capacity: (w(phi) - w(phi_prev)) /
            // (phi - phi_prev). With the lumped capacity matrix this conserves total
            // moisture across a step, unlike the tangent dw/dphi which creates/destroys
            // moisture where the sorption curve is steep (near saturation). See
            // SorptionSecantCapacity and doc/Moisture Governing Equations.md (D6).
            SorptionSecantCapacity secantCapacity(m_Material.sorptionCurve(), Variable::humidity);
            Cap(secantCapacity);
        }

        //////////////////////////////////////////////////////////////////////
        /// Functions for flux calculations
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeWaterLiquidTransportation()
           && m_Material.hasLiquidTransportationCurve())
        {
            CondFlux(LiquidTransportationCurve(m_Material.liquidTransportationCurve(), m_Material));
        }
    }

}   // namespace HygroThermFEM
